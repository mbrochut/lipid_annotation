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
        validation_dictionnary = {'name': 'PE (18:0/20:4)', 'abbreviation': 'PE', 'cls': None, 'sub_cls': None, 'prefix': '',
         'fa1': '18', 'u1': '0', 'fa2': '20', 'u2': '4', 'fa3': None, 'u3': None, 'add': '', 'match': True, 'known_chain': ''}
        self.assertDictEqual(validation_dictionnary,parsed_info)
        
    def test_valid_lipid_without_parentheses(self):
        l1 = Lipid(self.valid_lipid_without_parentheses)
        parsed_info = l1.parse_lipid_name()
        self.assertIsNotNone(parsed_info)
        validation_dictionnary = {'name': 'PE d18:0/20:4', 'abbreviation': 'PE', 'cls': None, 'sub_cls': None, 'prefix': 'd',
         'fa1': '18', 'u1': '0', 'fa2': '20', 'u2': '4', 'fa3': None, 'u3': None, 'add': '', 'match': True, 'known_chain': ''}
        self.assertDictEqual(validation_dictionnary,parsed_info)
        
    def test_invalid_lipid(self):
        l1 = Lipid(self.invalid_lipid)
        parsed_info = l1.parse_lipid_name()
        self.assertFalse(parsed_info['match'])

    

if __name__ == '__main__':
    unittest.main()