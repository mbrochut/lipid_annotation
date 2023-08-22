import unittest
from lipid_name import Lipid


class TestLipidParsing(unittest.TestCase):
    def setUp(self):
        self.valid_lipid_with_parentheses = "PE (18:0/20:4)"
        self.valid_lipid_without_parentheses = "PE d18:0/20:4"
        self.invalid_lipid = "InvalidLipidName"

    def test_valid_lipid_with_parentheses(self):
        l1 = Lipid(self.valid_lipid_with_parentheses)
        parsed_info = l1.parse_lipid_name()
        self.assertIsNotNone(parsed_info)
        validation_dictionnary = {'entry_name': 'PE', 'prefix': '',
         'fa1': '18', 'u1': '0', 'knownFA': '', 'fa2': '20', 'u2': '4', 'fa3': None, 'u3': None, 'additional_info': '', 'match': True}
        self.assertDictEqual(validation_dictionnary,parsed_info)
        
    def test_valid_lipid_without_parentheses(self):
        l1 = Lipid(self.valid_lipid_without_parentheses)
        parsed_info = l1.parse_lipid_name()
        self.assertIsNotNone(parsed_info)
        validation_dictionnary = {'entry_name': 'PE', 'prefix': 'd',
         'fa1': '18', 'u1': '0', 'knownFA': '', 'fa2': '20', 'u2': '4', 'fa3': None, 'u3': None, 'additional_info': '', 'match': True}
        self.assertDictEqual(validation_dictionnary,parsed_info)
        
    def test_invalid_lipid(self):
        l1 = Lipid(self.invalid_lipid)
        parsed_info = l1.parse_lipid_name()
        self.assertFalse(parsed_info['match'])

class TestLipidAttribut(unittest.TestCase):

    def test_attributs_with_valide_lipid(self):
        l1 = Lipid("LPE d18:0/20:4")
        lipid_attributs = l1.to_dict()
        validation_dictionnary = {'name': 'LPE d18:0/20:4', 'full_name': None, 'entry_name': 'LPE', 'cls': 'PE', 'sub_cls': "LPE",
         'fatty_acid_chains': ['18', '20', None], 'unsaturation': ['0', '4', None], 'knownFA': '', 'known_fatty_acid_chain': [], 'known_unsaturation': [],
          'prefix': 'd', 'additional_info': ''}
        self.assertDictEqual(validation_dictionnary,lipid_attributs)

    def test_attributs_with_FA_Known(self):
        l1 = Lipid("TAG 50:5_FA14:2+H20")
        lipid_attributs = l1.to_dict()
        print(lipid_attributs)
        validation_dictionnary = {'name': 'TAG 50:5_FA14:2+H20', 'full_name': None, 'entry_name': 'TAG', 'cls': 'TG', 'sub_cls': None,
         'fatty_acid_chains': ['50', None, None], 'unsaturation': ['5', None, None], 'knownFA': "FA", 'known_fatty_acid_chain': ["14"], 'known_unsaturation': ["2"],
          'prefix': "", 'additional_info': '+H20'}
        self.assertDictEqual(validation_dictionnary,lipid_attributs)




if __name__ == '__main__':
    unittest.main()