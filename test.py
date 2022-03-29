import unittest, json
from Bio import SeqIO
from path_variables import *

class TestHistonedbTypeClassifier(unittest.TestCase):
    def setUp(self):
        from prediction.histonedb_classifier import HistonedbTypeClassifier
        # self.classifier = HistonedbTypeClassifier()
        with open(VARIANTS_JSON) as f:
            variant_json = json.loads(f.read())
            self.classification_tree = variant_json['tree']
        self.classification_tree.pop('Archaeal')
        self.classification_tree.pop('Viral')
        self.classifier = HistonedbTypeClassifier(classification_tree=self.classification_tree)

    def test_create_hmms(self):
        # self.classifier.create_hmms()
        self.classifier.create_hmms(seed_directory=os.path.join(DATA_DIRECTORY, "draft_seeds"))
        self.assertEqual(self.classifier.combined_hmm_file, os.path.join(COMBINED_HMM_DIRECTORY, 'types_combined.hmm'))

    def test_predict_type(self):
        # self.classifier.create_hmms()
        self.classifier.create_hmms(seed_directory=os.path.join(DATA_DIRECTORY, "draft_seeds"))
        res = self.classifier.predict(sequences=os.path.join(PREDICTION_DIRECTORY, "test"))
        self.classifier.save_prediction_info(file_name=os.path.join(PREDICTION_RESULTS_DIRECTORY, "res_types.csv"))

        with open(os.path.join(PREDICTION_DIRECTORY, "test")) as t:
            test_data = list(SeqIO.parse(t, "fasta"))

        self.assertIsInstance(self.classifier.prediction_info, list)
        self.assertIsInstance(res, list)
        self.assertIsInstance(self.classifier.predicted_results, list)
        self.assertGreater(len(self.classifier.prediction_info), len(test_data))
        try: #NP_001295191.1
            self.assertEqual(len(res), len(test_data)-1)
            self.assertEqual(len(self.classifier.predicted_results), len(test_data) - 1)
        except:
            print(f"Symetric_difference is {set(res.values('accession')).symmetric_difference(set([s.id for s in test_data]))}")
            raise

    def test_dump_results(self):
        import pickle

        # self.classifier.create_hmms()
        # self.classifier.create_blastdbs()
        self.classifier.create_hmms(seed_directory=os.path.join(DATA_DIRECTORY, "draft_seeds"))
        res = self.classifier.predict(sequences=os.path.join(PREDICTION_DIRECTORY, "test"))
        self.classifier.dump_results(file_name=os.path.join(PREDICTION_RESULTS_DIRECTORY, "res_types.pickle"))

        with open(os.path.join(PREDICTION_RESULTS_DIRECTORY, "res_types.pickle"), 'rb') as f:
            data_new = pickle.load(f)

        self.assertEqual(len(data_new), len(self.classifier.prediction_info))
        self.assertEqual(len(data_new.get_keys()), len(self.classifier.prediction_info.get_keys()))

class TestHistonedbVariantClassifier(unittest.TestCase):
    def setUp(self):
        from prediction.histonedb_classifier import HistonedbVariantClassifier
        # self.classifier = HistonedbVariantClassifier()
        with open(VARIANTS_JSON) as f:
            variant_json = json.loads(f.read())
            self.classification_tree = variant_json['tree']
        self.classification_tree.pop('Archaeal')
        self.classification_tree.pop('Viral')
        self.classifier = HistonedbVariantClassifier(classification_tree=self.classification_tree)

    def test_create_hmms(self):
        # self.classifier.create_hmms()
        self.classifier.create_hmms(seed_directory=os.path.join(DATA_DIRECTORY, "draft_seeds"))
        self.assertEqual(self.classifier.combined_hmm_file, os.path.join(COMBINED_HMM_DIRECTORY, 'types_combined.hmm'))

    def test_create_blastdbs(self):
        self.classifier.create_blastdbs()

    def test_predict_variant(self):
        # self.classifier.create_hmms()
        self.classifier.create_hmms(seed_directory=os.path.join(DATA_DIRECTORY, "draft_seeds"))
        self.classifier.create_blastdbs()
        self.classifier.create_hmms(seed_directory=os.path.join(DATA_DIRECTORY, "draft_seeds"))
        res = self.classifier.predict(sequences=os.path.join(PREDICTION_DIRECTORY, "test"))
        self.classifier.save_prediction_info(file_name=os.path.join(PREDICTION_RESULTS_DIRECTORY, "res_variants.csv"))

        with open(os.path.join(PREDICTION_DIRECTORY, "test")) as t:
            test_data = list(SeqIO.parse(t, "fasta"))

        self.assertIsInstance(self.classifier.prediction_info, list)
        self.assertIsInstance(res, list)
        self.assertIsInstance(self.classifier.predicted_results, list)
        self.assertGreater(len(self.classifier.prediction_info), len(test_data))
        self.assertGreater(len(self.classifier.prediction_info),
                           len(list(filter(lambda d: d['best'], self.classifier.prediction_info))))
        try: #NP_001295191.1
            self.assertEqual(len(list(filter(lambda d: d['best'], self.classifier.prediction_info))), len(test_data)-1)
            self.assertEqual(len(res), len(test_data)-1)
            self.assertEqual(len(self.classifier.predicted_results), len(test_data)-1)
        except:
            print(f"Symetric_difference is {set(res.values('accession')).symmetric_difference(set([s.id for s in test_data]))}")
            raise

    def test_dump_results(self):
        import pickle

        # self.classifier.create_hmms()
        self.classifier.create_hmms(seed_directory=os.path.join(DATA_DIRECTORY, "draft_seeds"))
        self.classifier.create_blastdbs()
        res = self.classifier.predict(sequences=os.path.join(PREDICTION_DIRECTORY, "test"))
        self.classifier.dump_results(file_name=os.path.join(PREDICTION_RESULTS_DIRECTORY, "res_variants.pickle"))

        with open(os.path.join(PREDICTION_RESULTS_DIRECTORY, "res_variants.pickle"), 'rb') as f:
            data_new = pickle.load(f)

        self.assertEqual(len(data_new), len(self.classifier.prediction_info))
        self.assertEqual(len(data_new.get_keys()), len(self.classifier.prediction_info.get_keys()))

if __name__ == '__main__':
    unittest.main()