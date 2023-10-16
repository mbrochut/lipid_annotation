import pandas as pd
import requests

df_lipids = pd.read_excel("./results/lipid_t-test_with_lipid_name_convention.xlsx",index_col=0)


lipid_test = requests.get("https://www.swisslipids.org/api/search?term=Phosphatidate(36:2)")

print(lipid_test.json())

print(df_lipids)

def get_swissLipid_match(lipid_name):
    header = {"User-Agent": "Mozilla/5.0 (Windows; U; Windows NT 5.1; en-US; rv:1.9.2.8) Gecko/20100722 Firefox/3.6.8 GTB7.1 (.NET CLR 3.5.30729)", "Referer": "http://example.com"}
    proxy = {"http": "209.38.233.119:8080"}
    URL = "https://www.swisslipids.org/api/search?term={}".format(lipid_name)
    response = requests.get(URL,proxies=proxy, headers = header,timeout=30)
    swiss_id = None
    if response.status_code == 200:
        json_format = response.json()
        if len(json_format) > 1:
            print(json_format)
            for index, n in enumerate(json_format):
                if n['classification_level'] == "Class":
                    swiss_id = response.json()[index]['entity_id']
        else:
            swiss_id = response.json()[0]['entity_id']
        print(swiss_id)
        return swiss_id
    elif response.status_code >= 400:
        print("no response")
        return None
    

df_lipids['lipid_swiss_convention'] = df_lipids['lipid_name_convention'].str.replace(r"(^FFA)","FA")

df_lipids_day0 = df_lipids[df_lipids['day']=='D0']

df_lipids_day0_1 = df_lipids_day0.iloc[:20,:]

df_lipids_day0_1['SwissLipids_ID'] = df_lipids_day0_1['lipid_swiss_convention'].apply(lambda x: get_swissLipid_match(x))
print(df_lipids_day0_1)