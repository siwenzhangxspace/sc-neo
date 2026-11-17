#!/usr/bin/env python3
"""Regression tests for the scoring and evidence contracts; standard library only."""
import math
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from candidate_utils import peptide_hla_key
from candidate_expression_support import attach_expression_support
from integrate_candidate_results import index_iedb
from normalize_expression import normalize_rows, import_tpm
from score_candidates_topsis import rank_cohort, pci_module, RNA_REFERENCE
from scneo_utils import read_tsv, write_tsv


def example(**changes):
    row = dict(patient_id="P01", source_type="mutation", gene="TP53",
               peptide="SIINFEKLA", peptide_length="9", hla_allele="HLA-A*02:01",
               var_key="chr1_10_A_T", HGVSp_Short="p.A10V",
               mean_malignant_TPM=str(RNA_REFERENCE), expression_unit="TPM",
               PCI_DB_enabled="true", search_completed="true", pci_status="no_hit",
               healthy_blast_match="false", tumor_full_length_exact="false",
               max_TABR_BERT_rank="0.5", IEDB_immunogenicity_score="0",
               **{"Aff.nM": "150", "%Rank_EL": "0.5"})
    row.update(changes)
    return row


class ScoringTests(unittest.TestCase):
    def test_fixed_ideals_have_known_answer(self):
        # P=C=T=I=.5, R=1; distances calculated independently.
        row = rank_cohort([example()])[0]
        d0 = math.sqrt(.125**2 + .1**2 + .125**2 + .15**2 + .075**2)
        d1 = math.sqrt(.125**2 + .1**2 + .125**2 + .075**2)
        self.assertAlmostEqual(float(row["TOPSIS_score"]), 100*d0/(d0+d1), places=12)

    def test_four_modules_without_tcr(self):
        row = rank_cohort([example(max_TABR_BERT_rank="")])[0]
        d0 = math.sqrt(.2**2 + .15**2 + .15**2 + .075**2)
        d1 = math.sqrt(.2**2 + .15**2 + .075**2)
        self.assertEqual(row["weight_policy"], "TCR_missing_40_30_15_15")
        self.assertAlmostEqual(float(row["TOPSIS_score"]), 100*d0/(d0+d1), places=12)

    def test_missing_both_not_scored(self):
        row = rank_cohort([example(PCI_DB_enabled="false", max_TABR_BERT_rank="")])[0]
        self.assertEqual(row["TOPSIS_score"], "")
        self.assertEqual(row["top_decile"], "false")

    def test_pci_disabled_with_tcr_is_explicit_mode(self):
        row = rank_cohort([example(PCI_DB_enabled="false")])[0]
        self.assertEqual(row["weight_policy"], "PCI_missing_35_35_15_15")
        self.assertEqual(row["score_eligible"], "true")

    def test_zero_expression_is_measured_not_missing(self):
        zero, missing = [rank_cohort([example(mean_malignant_TPM=x)])[0] for x in ("0","")]
        self.assertEqual(zero["topsis_rna_expression"], "0")
        self.assertNotEqual(zero["TOPSIS_score"], "")
        self.assertEqual(missing["TOPSIS_score"], "")

    def test_scores_independent_of_cohort(self):
        one = rank_cohort([example()])[0]["TOPSIS_score"]
        many = rank_cohort([example(), example(gene="EGFR", mean_malignant_TPM="0")])
        self.assertEqual(one, many[0]["TOPSIS_score"])

    def test_top_decile_budget_and_stable_ties(self):
        rows = [example(gene=f"G{i:02}") for i in range(11)]
        a = rank_cohort(rows)
        b = rank_cohort(list(reversed(rows)))
        self.assertEqual({r["gene"] for r in a if r["top_decile"]=="true"}, {"G00","G01"})
        self.assertEqual({r["gene"]:r["cohort_rank"] for r in a},
                         {r["gene"]:r["cohort_rank"] for r in b})

    def test_duplicate_identity_fails(self):
        with self.assertRaises(ValueError):
            rank_cohort([example(),example()])

    def test_no_hit_is_neutral_and_short_hit_is_attenuated(self):
        self.assertEqual(pci_module(example()), .5)
        row = example(pci_status="hit", blast_identity="100", blast_evalue="0", query_coverage=".5")
        self.assertEqual(pci_module(row), .75)
        self.assertEqual(pci_module(dict(row, tumor_full_length_exact="true")), 1)

    def test_failed_search_is_not_no_hit(self):
        with self.assertRaises(ValueError):
            rank_cohort([example(search_completed="false")])

    def test_legacy_expression_requires_regeneration(self):
        row = example(RNA_expression_support="5")
        del row["mean_malignant_TPM"]
        with self.assertRaises(ValueError):
            rank_cohort([row])

    def test_invalid_rank_and_incomplete_prediction_fail(self):
        for updates in ({"max_TABR_BERT_rank":"1.1"}, {"tabr_bert_status":"missing_score"},
                        {"mean_malignant_TPM":"NaN"}, {"healthy_blast_match":"true"}):
            with self.subTest(updates=updates), self.assertRaises(ValueError):
                rank_cohort([example(**updates)])

    def test_hla_aware_join_preserves_different_scores(self):
        rows = [dict(peptide="SIINFEKLA", hla_allele=hla, immunogenicity_score=score)
                for hla,score in [("HLA-A*02:01",".1"),("HLA-A*24:02",".2")]]
        scores = index_iedb(rows)
        self.assertEqual(scores[peptide_hla_key(dict(peptide="SIINFEKLA",hla_allele="HLA-A2402"))]["immunogenicity_score"], "0.2")
        self.assertEqual(len(scores),2)

    def test_mean_includes_all_malignant_cells_and_zeros(self):
        expr = [dict(patient_id="P01",cell_id=c,gene_name="TP53",tpm=v)
                for c,v in [("C1","100"),("C2","0")]]
        manifest = [dict(patient_id="P01",cell_id=c,cell_type="tumor") for c in ["C1","C2"]]
        row = attach_expression_support([example(supporting_tumor_cell_ids="C1")],expr,manifest)[0]
        self.assertEqual(float(row["mean_malignant_TPM"]), 50)
        self.assertEqual(row["expression_detected_cells"],1)
        missing = attach_expression_support([example()],expr[:1],manifest)[0]
        self.assertEqual(missing["mean_malignant_TPM"],"")
        self.assertEqual(missing["expression_mapping_status"],"incomplete_gene_measurement")

    def test_gene_length_normalization_and_empty_library(self):
        rows = [dict(patient_id="P01",cell_id="C1",gene_id=g,gene_name=g,raw_count="10",gene_length=length)
                for g,length in [("G1","1000"),("G2","2000")]]
        normalized = normalize_rows(rows)
        self.assertAlmostEqual(float(normalized[0]["tpm"]), 2e6/3)
        self.assertAlmostEqual(sum(float(r["tpm"]) for r in normalized),1e6)
        empty = normalize_rows([dict(r,raw_count="0") for r in rows])
        self.assertTrue(all(r["tpm"]=="" for r in empty))

    def test_import_rejects_negative_and_duplicate_measurements(self):
        row = dict(patient_id="P01",cell_id="C1",gene_name="TP53",tpm="0")
        self.assertEqual(import_tpm([row])[0]["tpm"],"0")
        for invalid in ([dict(row,tpm="-1")],[row,row]):
            with self.assertRaises(ValueError):
                import_tpm(invalid)

    def test_tpm_import_and_candidate_support_cli(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp=Path(directory)
            rows=[dict(patient_id='P01',cell_id=c,gene_name='TP53',tpm=tpm) for c,tpm in [('C1','100'),('C2','0')]]
            write_tsv(tmp/'tpm.tsv',rows,list(rows[0]))
            self.cli('normalize_expression.py','--tpm',tmp/'tpm.tsv','--out-long',tmp/'long.tsv',
                     '--out-raw-matrix',tmp/'raw.tsv','--out-normalized-matrix',tmp/'norm.tsv','--out-tpm-matrix',tmp/'matrix.tsv')
            write_tsv(tmp/'pmhc.tsv',[example()],list(example()))
            manifest=[dict(patient_id='P01',cell_id=c,cell_type='tumor') for c in ['C1','C2']]
            write_tsv(tmp/'manifest.tsv',manifest,list(manifest[0]))
            self.cli('candidate_expression_support.py','--pmhc',tmp/'pmhc.tsv','--expression',tmp/'long.tsv',
                     '--manifest',tmp/'manifest.tsv','--out',tmp/'support.tsv')
            self.assertEqual(float(read_tsv(tmp/'support.tsv')[0]['mean_malignant_TPM']),50)

    def test_no_tcr_cli_does_not_call_predictor(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp=Path(directory)
            write_tsv(tmp/'pmhc.tsv',[example()],list(example()))
            self.cli('predict_tcr_recognition.py','--pmhc',tmp/'pmhc.tsv','--cmd','nonexistent_predictor', '--out',tmp/'tcr.tsv')
            row=read_tsv(tmp/'tcr.tsv')[0]
            self.assertEqual(row['tabr_bert_status'],'no_valid_tcr')
            self.assertEqual(row['max_TABR_BERT_rank'],'')

    def test_cd4_rows_do_not_enter_tcr_prediction(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp=Path(directory)
            write_tsv(tmp/'pmhc.tsv',[example()],list(example()))
            row=dict(patient_id='P01',cell_id='T1',cdr3_beta_aa='CASSLGQETQYF',trust4_status='ok',tcell_subtype='CD4')
            write_tsv(tmp/'tcr.tsv',[row],list(row))
            self.cli('predict_tcr_recognition.py','--pmhc',tmp/'pmhc.tsv','--tcr',tmp/'tcr.tsv',
                     '--cmd','nonexistent_predictor','--out',tmp/'out.tsv')
            self.assertEqual(read_tsv(tmp/'out.tsv')[0]['tabr_bert_status'],'no_valid_tcr')

    def cli(self,script,*args):
        result=subprocess.run([sys.executable,'-B',str(Path(__file__).parent/script),*map(str,args)],capture_output=True,text=True)
        self.assertEqual(result.returncode,0,result.stdout+'\n'+result.stderr)


if __name__ == "__main__":
    unittest.main()
