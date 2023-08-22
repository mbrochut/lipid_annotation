import re
import pandas as pd
import warnings

class Lipid:
    def __init__(self, name):
        self.name = name
        self.full_name = None # to convert lipid class to full name ex: PE = Phosphatidylethanolamine
        self.entry_name = None
        self.cls = None
        self.sub_cls = None
        self.fatty_acid_chains = []
        self.unsaturation = []
        self.known_fatty_acid_chain = []
        self.known_unsaturation = []
        self.prefix = None
        self.additional_info = None

    # set and add functions
    def add_fatty_acid_chain(self, chain, unsaturation):
        if len(self.fatty_acid_chains) < 3:
            self.fatty_acid_chains.append(chain)
            self.unsaturation.append(unsaturation)
        else:
            print("Maximum number of fatty acid chains (3) reached.")

    def add_known_fatty_acid_chain(self, chain, unsaturation):
        self.known_fatty_acid_chain.append(chain)
        self.known_unsaturation.append(unsaturation)
  
    
    def from_name_to_attributes(self):
        parsed_name = self.parse_lipid_name()

        #add carbon chain
        for index, n in enumerate([1,2,3]):
            self.add_fatty_acid_chain(parsed_name[f'fa{n}'],parsed_name[f'u{n}'])
        
        self.update_attributes(parsed_name)
        self.check_lipide_class()
        self.check_lipid_sub_clas()

    def parse_lipid_name(self):
        lipid_pattern = re.compile(
            r"^(?P<entry_name>[A-Za-z_]+)[\(\s]*"  # Abbreviation (letters and underscores), followed by optional open parenthesis or spaces
            r"(?P<prefix>[mdte]*|[-PO]*)"  # Prefix (optional characters), can be 'm', 'd', 't', 'e', '-', 'P', or 'O'
            r"(?P<fa1>[0-9]+):(?P<u1>[0-9]+)[/|_]?(?P<known_chain>[FA]*)?"  # Chain 1: Carbon chain length and unsaturation count, followed by optional slash or underscore and known chain designation
            r"((?P<fa2>[0-9]+):(?P<u2>[0-9]+))*[/|_]?((?P<fa3>[0-9]+):(?P<u3>[0-9]+))*[\)\s]*"  # Optional Chain 2 and Chain 3 (similar structure as Chain 1), enclosed in parentheses
            r"(?P<add>.*)"  # Additional information (anything remaining)
        )
        result = {
        'name': self.name,
        'entry_name': None,
        'cls': None,
        'sub_cls': None,
        'prefix': None,
        'fa1': None,
        'u1': None,
        'fa2': None,
        'u2': None,
        'fa3': None,
        'u3': None,
        'add': None,
        'match': None
        }

        # parse the name using regex
        l_res = lipid_pattern.match(self.name)
        if l_res == None:
            warnings.warn(self.name,SyntaxWarning)
            regex_dict = {"match" : False}
        else:
            regex_dict = l_res.groupdict()
            regex_dict['match'] = True
        
        result.update(regex_dict)
        return result

    def update_attributes(self, update_dict):
        for key, value in update_dict.items():
            setattr(self, key, value)

    def check_lipide_class(self):
        string_to_class = {
            r'(ce|lce)+$': 'CE',
            r'(pe|lpe)+$': 'PE',
            r'(hcer|hexcer)+$': 'HexCer',
            r'(cer|lcer)+$': 'Cer',
            r'(MAG|MG)+$': "MG",
            r'(DAG|DG)+$': "DG",
            r'(TAG|TG)+$': "TG",
            r'(SM)+$': "SM",
            r'(DCER)+$': "DhCer",
            r'(pc|lpc)+$': "PC",
            r'(pi|lpi)+$': "PI",
            r'(pg)+$': "PG",
            r'(ps)+$': "PS",
            r'.*(ic)+$': "FFA" #if name end in "ic" it's a Free Fatty acid


            # Add more mappings here as needed
        }

        string_class = self.entry_name
        if string_class == None:
            return
        
        items_per_indicies =[items for key, items in string_to_class.items()]
        dictkeys_pattern = re.compile('|'.join(string_to_class), re.IGNORECASE)
        match = re.match(dictkeys_pattern, string_class)
        
        if match == None:
            print("class is unkown for",self.name)
            return
        
        index = match.lastindex - 1
        l_class = items_per_indicies[index]
        self.cls = l_class

    def check_lipid_sub_clas(self):
        #If the lipid is not a free fatty acid (FFA) it will check if the name begin with L for liso Lipids.
        if self.cls != "FFA" and self.entry_name.lower().startswith(r"l"):
            self.sub_cls = "L"+self.cls

    
    
    def display_info(self):
        print(f"Name: {self.name}")
        print(f"Class: {self.cls}")
        print(f"Subclass: {self.sub_cls}")
        print("Fatty Acid Chains:")
        for idx, chain in enumerate(self.fatty_acid_chains):
            print(f"  Chain {idx + 1}: {chain} - Unsaturation: {self.unsaturation[idx]}")
        print(f"Prefix: {self.prefix}")
        print(f"Additional Info: {self.additional_info}")

    def short_info(self):
        print(f"Name: {self.name}")


def check_known_FA_in_name(matches):
    if matches['known_chain'] == 'FA':
        matches['known_chain'] == True
        matches['fa1'] = matches['fa2']
        matches['kfa1'] = matches['u2']
        matches["fa2"] = None
        matches["u2"] = None


def check_lipide_class(matches):
    string_to_class = {
        r'(ce|lce)+$': 'CE',
        r'(pe|lpe)+$': 'PE',
        r'(hcer|hexcer)+$': 'HexCer',
        r'(cer|lcer)+$': 'Cer',
        r'(MAG|MG)+$': "MG",
        r'(DAG|DG)+$': "DG",
        r'(TAG|TG)+$': "TG",
        r'(SM)+$': "SM",
        r'(DCER)+$': "DhCer",
        r'(pc|lpc)+$': "PC",
        r'(pi|lpi)+$': "PI",
        r'(pg)+$': "PG",
        r'(ps)+$': "PS",
        r'.*(ic)+$': "FFA" #if name end in "ic" it's a Free Fatty acid


        # Add more mappings here as needed
    }

    string_class = matches['abbreviation']
    if string_class == None:
        return
    
    items_per_indicies =[items for key, items in string_to_class.items()]
    dictkeys_pattern = re.compile('|'.join(string_to_class), re.IGNORECASE)
    match = re.match(dictkeys_pattern, string_class)
    
    if match == None:
        print("class is unkown for",matches['name'])
        matches["class"] = "unknown"
        return
    index = match.lastindex - 1
    l_class = items_per_indicies[index]
    matches["class"] = l_class
    if l_class != "FFA":
        check_lipid_sub_clas(matches)
    
def check_lipid_sub_clas(matches):
    if matches['name_class'].lower().startswith(r"l"):
        matches['sub-class'] = "L"+matches['class']

def convention_database(name_DB,cls,sub_cls,prefix,sep,parenthesis,kfc,**kwargs):
    if not name_DB in ['swiss','CHUV']:
        return cls,prefix,sep,parenthesis,kfc
    
    if name_DB == "swiss":
        if cls == "FFA":
            cls = "FA"
        elif cls in ["Cer","HexCer","SM","PC"]:
            prefix = "d"
            sep = "/"
        elif cls == "DhCer":
            cls = "Cer"
            prefix = "d"
            sep = "/"
        return cls,prefix,sep,parenthesis,kfc
    
    if name_DB == "CHUV":
        if cls in ["Cer","HexCer","DhCer"]:
            prefix = "d"
            sep = "/"
        parenthesis = False
        if cls == "TG":
            kfc = True
        if sub_cls in ["LPC","LPE"]:
            cls = sub_cls
        if cls == "FFA":
            print(kwargs)
            cls = kwargs['name']
        return cls,prefix,sep,parenthesis,kfc

def format_for_lipid_name(cls, prefix, fc1,fu1,fc2,fu2,fc3,fu3,sep = "_",parenthesis = True,sub_cls=None,name_DB="",kfc = False,**kwargs):
    
    #arguments = {"cls":cls,"prefix": prefix,"fc1":fc1,"fc2":fc2,"fc3":fc3,"fu1":fu1,"fu2":fu2,"fu3":fu3,"sep" : sep}
    cls,prefix,sep,parenthesis,kfc = convention_database(name_DB,cls,sub_cls,prefix,sep,parenthesis,kfc,**kwargs)
    #carbon_chain
    carbon_chain = fc1+":"+fu1
    if not fc2 == None:
        carbon_chain = carbon_chain +sep +fc2+":"+fu2
    if not fc3 == None:
        carbon_chain = carbon_chain +sep +fc3+":"+fu3

    if kfc:
        carbon_chain = carbon_chain +sep+ kwargs['kfc1']+":"+kwargs['kfu1']

    if not sub_cls == None:
        cls = sub_cls
    if parenthesis:
        lipid_name = cls +"("+ prefix + carbon_chain +")"
    else:
        lipid_name = cls +" "+ prefix + carbon_chain
    return lipid_name


if __name__ == "__main__":

    # Example usage
    lipid = Lipid("PE (18:0/20:4)")
    lipid.cls = "Phosphatidylethanolamine"
    lipid.add_fatty_acid_chain("18:0", 0)
    lipid.add_fatty_acid_chain("20:4", 4)
    lipid.display_info()
    lipid.short_info()

    
    res = lipid.parse_lipid_name()
    print(res)

    lipid2 = Lipid("wfweoi (d18:0/20:4)")
    res2 = lipid2.parse_lipid_name()
    lipid2.from_name_to_attributes()
    lipid2.display_info()

    