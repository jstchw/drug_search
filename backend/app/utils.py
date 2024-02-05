import json
from Bio import Entrez
import re
import time


def format_json_drug(drug, product_name=None):
    calculated_properties = drug.get('calculated-properties', {}).get('property', [])

    groups = drug.get('groups', {}).get('group', [])
    groups = [groups] if isinstance(groups, str) else groups

    return {
        'name': drug.get('name'),
        'half_life': drug.get('half-life'),
        'classification': drug.get('classification', {}).get('class'),
        'groups': groups,
        'brands': [brand.get('name') for brand in drug.get('international-brands', {}).get('international-brand', [])] if isinstance(drug.get('international-brands', {}).get('international-brand'), list) else None,
        'iupac': next((prop['value'] for prop in calculated_properties if prop.get('kind') == 'IUPAC Name'), None),
        'formula': next((prop['value'] for prop in calculated_properties if prop.get('kind') == 'Molecular Formula'), None),
        'indication': drug.get('indication'),
        'product': product_name,
    }


def search_json(params, json_file_path, limit=10):
    start_time = time.time()
    matched_entries = []
    search_fields = []

    if params['search_type'] in ['generic_name', 'brand_name']:
        search_fields = ['drug']  # Search in 'drug' for both generic and brand names
    elif params['search_type'] == 'side_effect':
        search_fields = ['effect', 'treatment_disorder']  # Search in both for side effects

    gender_specified = params.get('sex') is not None and params['sex'].strip() != ''
    gender_filter = params.get('sex').lower() if gender_specified else None

    with open(json_file_path, 'r') as file:
        data = json.load(file)  # Load the entire JSON array

    for entry in data:
        # Filter by gender
        if gender_specified:
            if not any(gender_filter == g.lower() for g in entry.get('gender', [])):
                continue

        # Relaxed search:
        # Match even if more drugs are present in the entry
        if all(any(term.lower() in item.lower() for field in search_fields for item in entry.get(field, [])) for term in params['terms']):
            matched_entries.append(entry)

    # Sort the entries by publication date in descending order and apply the limit
    sorted_entries = sorted(matched_entries, key=lambda x: x['pub_date'], reverse=True)[:limit]

    print(f"Search took {time.time() - start_time:.5f} seconds", flush=True)
    return sorted_entries


def get_pubmed_metadata(pubmed_ids):
    id_string = ','.join(pubmed_ids)

    Entrez.email = "k21002495@kcl.ac.uk"
    handle = Entrez.efetch(db="pubmed", id=id_string, retmode="xml")

    records = Entrez.read(handle, validate=True)
    records = records['PubmedArticle']
    handle.close()

    return filter_pubmed_metadata(records)


def filter_pubmed_metadata(records):
    filtered_records = []

    for record in records:
        title = record['MedlineCitation']['Article']['ArticleTitle']

        title = re.sub(r'<.*?>', '', title)  # Remove HTML tags

        abstract = record['MedlineCitation']['Article']['Abstract']['AbstractText'][0] \
            if 'Abstract' in record['MedlineCitation']['Article'] else ''

        abstract = re.sub(r'<.*?>', '', abstract)  # Remove HTML tags

        # Extract authors
        authors = []
        author_list = record.get('MedlineCitation', {}).get('Article', {}).get('AuthorList', [])
        if isinstance(author_list, list):
            for author in author_list:
                # Ensure 'author' is a dictionary before attempting to access keys
                if isinstance(author, dict):
                    fore_name = author.get('ForeName', '')
                    last_name = author.get('LastName', '')
                    authors.append(f"{fore_name} {last_name}".strip())

        article_url = record['MedlineCitation']['PMID']

        article_date = record['MedlineCitation']['Article']['ArticleDate']
        if article_date:
            article_date = article_date[0]
            article_date = f"{article_date['Year']}-{article_date['Month']}-{article_date['Day']}"
        else:
            article_date = ''

        country = record['MedlineCitation']['MedlineJournalInfo']['Country']

        filtered_records.append({
            'title': title,
            'abstract': abstract,
            'authors': authors,
            'date': article_date,
            'country': country,
            'url': f"https://pubmed.ncbi.nlm.nih.gov/{article_url}"
        })

    return filtered_records

