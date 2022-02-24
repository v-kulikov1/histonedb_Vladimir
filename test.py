import unittest

class TestStringMethods(unittest.TestCase):
    def setUp(self):
        from prediction.histonedb_classifier import HistonedbClassifier
        # self.classifier = HistonedbClassifier()
        self.classifier = HistonedbClassifier(classification_tree={'H2A':{'H2A.B': 'null', 'H2A.Z': 'null', 'macroH2A': 'null'}})

    def test_create_hmms(self):
      self.classifier.create_hmms()

    # def test_create_blastdbs(self):
    #   self.assertTrue('FOO'.isupper())
    #   self.assertFalse('Foo'.isupper())


if __name__ == '__main__':
    unittest.main()