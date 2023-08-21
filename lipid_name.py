import re
import pandas as pd


def pars_lipide_name(name):
    name = name.replace(" ","")
    l_pat = re.compile(
        r"^(?P<name_class>[A-Za-z123-_]+)\((?P<prefix>[mdte]*|[-PO]*)(?P<fc1>[0-9]+):(?P<fu1>[0-9]+)[/|_]*(?P<known_chain>[FA]*)((?P<fc2>[0-9]+):(?P<fu2>[0-9]+))*[/|_]*((?P<fc3>[0-9]+):(?P<fu3>[0-9]+))*\)(?P<add>.*)")
    result = {
    'name': name,
    'name_class': None,
    'class': None,
    'sub-class': None,
    'prefix': None,
    'fc1': None,
    'fu1': None,
    'fc2': None,
    'fu2': None,
    'fc3': None,
    'fu3': None,
    'add': None,
    'known_chain': None,
    'kfc1': None,
    'kfu1': None,
    'match': None
    }
    """The sphingoid backbone is annotated by the number of hydroxyl groups in the sphingoid base
    (m for mono, d for di, t for tri) and separated by a slash from the number of carbons:number 
    of double bonds of the N-linked fatty acid. Positions of hydroxyl groups and double bonds including
    geometry are indicated as described for fatty acyls (FA)."""

    # parse the name using regex
    # parse the name using regex
    l_res = l_pat.match(name)
    if l_res == None:
        print("impossible to parse name: ", name)
        regex_dict = {"match" : False}
    else:
        regex_dict = l_res.groupdict()
        regex_dict['match'] = True
    
    result.update(regex_dict)
    check_FA_in_name(result)
    check_lipide_class(result)

    
    return result

def check_FA_in_name(matches):
    if matches['known_chain'] == 'FA':
        matches['known_chain'] == True
        matches['kfc1'] = matches['fc2']
        matches['kfu1'] = matches['fu2']
        matches["fc2"] = None
        matches["fu2"] = None


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
        r'.*(ic)+$': "FFA"


        # Add more mappings here as needed
    }

    string_class = matches['name_class']
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

    print("Lipid Annotation Programme")
    df = pd.read_excel("./Data/lipid_name_testing.xlsx")
    print(df)

    