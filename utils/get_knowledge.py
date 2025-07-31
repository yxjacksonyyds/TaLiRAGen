from chembl_webresource_client.new_client import new_client
import requests

drug = new_client.drug
target = new_client.target
activity = new_client.activity

def get_id(pref_name):
    if target.filter(pref_name__iexact= pref_name).filter(organism='Homo sapiens'):
        uniprotID = target.filter(pref_name__iexact= pref_name).filter(organism='Homo sapiens')[0]['target_components'][0]['accession']
        chemblid = target.filter(pref_name__iexact= pref_name).filter(organism='Homo sapiens')[0]['target_chembl_id']
        return uniprotID,chemblid
    else:
        print("No target found")
        
#通过靶点信息获取uniprot信息

    # 提取信息的函数

def get_icnm(chemblid):
    herg_activities = activity.filter(target_chembl_id=chemblid).filter(standard_type="IC50").order_by(["standard_value"])
    list = []
    for i in herg_activities[:5]:
        # 调用函数
        list = extract_information(i,list)
    return list

def get_drug(chemid):
    herg_activities = drug.filter(target_chembl_id=chemid)
    list = []
    for i in herg_activities[:5]:
        # 调用函数
        list = extract_information(i,list)
    return list
    
def extract_information(entry,list):
    smiles = entry.get('canonical_smiles')
    ki_value = entry.get('standard_value')
    chembl_id = entry.get('assay_chembl_id')
    description = entry.get('assay_description')
    
    # print(f"分子SMILES: {smiles}")
    # print(f" Ki值: {ki_value} nM")
    # print(f" ChEMBL ID: {chembl_id}")
    # print(f"描述：{description}")
    # print()  # 打印空行以分隔不同的条目

    dict = {"分子SMILES":smiles,"Ki（越高说明结合能力越弱）":ki_value,"CHEMBLID":chembl_id,"description":description}
    list.append(dict)
    
    return list


def get_activity(chemblid):
    herg_activities = activity.filter(target_chembl_id=chemblid).filter(standard_type="Ki").order_by(["standard_value"])
    list = []
    for i in herg_activities[:5]:
        # 调用函数
        list = extract_information(i,list)
    return list

def get_bind(uniprotid):
    ip = f'https://rest.uniprot.org/uniprotkb/{uniprotid}.json?fields=ft_binding%2Cft_site%2Cft_act_site'
    response = requests.get(ip)
    if response.status_code == 200:
        data = response.json()
        return data
    else:
        print("Invalid uniprot ID")

def get_pref_name(chemblid):
    # 设置 URL
    url = f"https://www.ebi.ac.uk/chembl/api/data/target/{chemblid}.json"
    
    # 发送 GET 请求
    response = requests.get(url)
    
    # 检查请求是否成功
    if response.status_code == 200:
        # 解析 JSON 数据
        data = response.json()
        return data['pref_name']
    else:
        return None

def get_text(key_word,max_result=2):
    from pubmed_sdk import PubMed
    pubmed = PubMed()
    results = pubmed.search(term=key_word)
    id_list = results['id_list']    # ['33725716', '33725717']
    details = pubmed.fetch_details(id_list).get('PubmedArticle')
    abstract = []
    for detail in details[:2]:
        article = detail['MedlineCitation']['Article']
        abstract.append(article['Abstract']['AbstractText'])
    return abstract

def get_from_uniprot_name(uniprot_name):
    ip = f'http://www.uniprot.org/uniprot/{uniprot_name}.json'
    response = requests.get(ip)
    if response.status_code == 200:
        data = response.json()

        uniprotid = data['primaryAccession']
        chem = data['uniProtKBCrossReferences']
        chembl_entries = [entry for entry in chem if entry.get("database") == "ChEMBL"]
        if len(chembl_entries)!=0:
            chemblid = chembl_entries[0]['id']
        else:
            return None
        return uniprotid,chemblid
    else:
        print("Invalid uniprot ID")
        return None


if __name__ == '__main__':

    pref_name = input("Please give me the pref_name of target:")
    uniprotid,chemblid = get_id(pref_name)
    get_activity(chemblid)