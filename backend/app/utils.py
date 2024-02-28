import json
from Bio import Entrez
import re
import time
from constants import age_groups, sex_groups


def format_json_drug(drug, product_name=None):
    """
    Format the given drug data into a standardized JSON format.
    
    Args:
        drug: dict, the drug data to be formatted.
        product_name: str, optional, the name of the product.
        
    Returns:
        dict: the formatted drug data in a standardized JSON format.
    """
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


def strict_search(entry, params, search_fields):
    """
    Perform a strict search on the given entry using the provided parameters and search fields.
    Currently NOT IN USE.

    Args:
        entry (dict): The entry to be searched.
        params (dict): The parameters for the search.
        search_fields (list): The fields to be searched in the entry.

    Returns:
        list: A list of matched entries.
    """
    matched_entries = []

    if all(re.search(r'(?:\b|\s){}(?:\b|\s|$)'.format(re.escape(term.lower())), item.lower()) for field in
            search_fields for item in entry.get(field, []) for term in params['terms']):
        matched_entries.append(entry)

    return matched_entries


def search_json(params, data, limit=10) -> list:
    """
    Search for entries in a JSON file based on the given parameters.

    :param params: A dictionary containing search parameters such as 'search_type', 'sex', 'age', 'country', 'search_mode', and 'terms'.
    :param json_file_path: The file path to the JSON file to be searched.
    :param limit: An integer specifying the maximum number of matched entries to return (default is 10).
    :return: A list of matched entries sorted by publication date in descending order.
    """
    start_time = time.time()
    matched_entries = []
    search_fields = []

    if params['search_type'] in ['generic_name', 'brand_name']:
        search_fields = ['drug']  # Search in 'drug' for both generic and brand names
    elif params['search_type'] == 'side_effect':
        search_fields = ['effect', 'treatment_disorder']  # Search in both for side effects

    # Check if gender actually exists in the parameters and is not empty
    gender_specified = params.get('sex') is not None and params.get('sex').strip() != ''
    gender_filter = params.get('sex').lower() if gender_specified else None

    # Check if age actually exists in the parameters and is not empty
    age_param = params.get('age')
    if isinstance(age_param, dict):
        age_specified = any(value is not None for value in age_param.values())
    else:
        age_specified = False

    country_specified = params.get('country') is not None and params.get('country').strip() != ''
    country_filter = params.get('country').lower() if country_specified else None

    for entry in data:
        # Gender filtering section -------------------------------------------
        if gender_specified:
            entry_gender = entry.get('gender')

            if entry_gender is None or gender_filter != entry_gender:
                continue

        # Age filtering section ----------------------------------------------
        if age_specified:
            entry_age = entry.get('age')

            # If the age is not specified in the entry, skip it
            if entry_age is None:
                continue

            # Init flags for min and max age
            matches_min_age = False
            matches_max_age = False

            # Checks for minimum age param
            if age_param.get('min') is not None:
                min_age = int(age_param['min'])

                # For the case where age is a list (range)
                if isinstance(entry_age, list):
                    # Check if the minimum age is less than or equal to the maximum age in the dataset
                    matches_min_age = min_age <= entry_age[1]
                else:
                    # Check if the minimum age is less than or equal to the age in the dataset (single value)
                    matches_min_age = min_age <= entry_age

            # Checks for maximum age param
            if age_param.get('max') is not None:
                max_age = int(age_param['max'])

                # For the case where age is a list (range)
                if isinstance(entry_age, list):
                    # Check if the maximum age is greater than or equal to the minimum age in the dataset
                    matches_max_age = max_age >= entry_age[0]
                else:
                    # Check if the maximum age is greater than or equal to the age in the dataset (single value)
                    matches_max_age = max_age >= entry_age

            # If the entry does not match the age range, skip it
            if (age_param.get('min') is not None and not matches_min_age) or \
                    (age_param.get('max') is not None and not matches_max_age):
                continue

        # Country filtering section ------------------------------------------
        if country_specified:
            entry_country = entry.get('country')

            if entry_country is None or country_filter != entry_country.lower():
                continue




        if params['search_mode'] == 'strict':
            # Strict search:
            # Match only if all drugs or side effects are present in the entry
            # NOW WORKS ONLY WITH GENERIC NAME OR BRAND NAME
            # DOES NOT WORK WITH SIDE EFFECTS!!!
            for field in search_fields:
                items = entry.get(field, [])
                if all(re.search(r'(?:\b|\s){}(?:\b|\s|$)'.format(re.escape(term.lower())), item.lower()) for item in
                       items for term in params['terms']):
                    # Check if the drug is exclusive in the entry
                    if len(items) == len(params['terms']):
                        matched_entries.append(entry)
        else:
            # Relaxed search:
            # Match even if more drugs or side effects are present in the entry
            # Regex explanation:
            # (?:\b|\s) - Match a word boundary or a whitespace character
            # {} - The search term
            # (?:\b|\s|$) - Match a word boundary, a whitespace character or the end of the string
            if all(any(re.search(r'(?:\b|\s){}(?:\b|\s|$)'.format(re.escape(term.lower())), item.lower()) for field in
                    search_fields for item in entry.get(field, [])) for term in params['terms']):
                matched_entries.append(entry)


    # Sort the entries by publication date in descending order and apply the limit
    sorted_entries = sorted(matched_entries, key=lambda x: x['pub_date'], reverse=True)[:limit]

    print(f"Search took {time.time() - start_time:.5f} seconds", flush=True)
    return sorted_entries


def get_pubmed_metadata(pubmed_ids):
    """
    Retrieves metadata from PubMed for the given PubMed IDs.

    Args:
        pubmed_ids (list): A list of PubMed IDs.

    Returns:
        list: A filtered list of PubMed metadata.
    """
    id_string = ','.join(pubmed_ids)

    Entrez.email = "k21002495@kcl.ac.uk"
    handle = Entrez.efetch(db="pubmed", id=id_string, retmode="xml")

    records = Entrez.read(handle, validate=True)

    if records is not None:
        records = records['PubmedArticle']
    handle.close()

    return filter_pubmed_metadata(records)


def filter_pubmed_metadata(records):
    """
    Generate a list of filtered PubMed metadata records based on the input records.
    Take only what's needed for the frontend.

    Args:
        records (list): A list of PubMed metadata records.

    Returns:
        list: A list of filtered PubMed metadata records containing title, abstract, authors, date, country, and URL.
    """
    filtered_records = []

    for record in records:
        title = record['MedlineCitation']['Article']['ArticleTitle']

        title = re.sub(r'<.*?>', '', title)  # Remove HTML tags

        abstract_section = record['MedlineCitation']['Article']['Abstract']['AbstractText'] \
            if 'Abstract' in record['MedlineCitation']['Article'] and \
            'AbstractText' in record['MedlineCitation']['Article']['Abstract'] and \
            len(record['MedlineCitation']['Article']['Abstract']['AbstractText']) > 0 else ''
        
        abstract = ''

        # Check if abstract_section is directly a string
        if isinstance(abstract_section, str):
            abstract = abstract_section
        # Check if abstract_section is a list
        elif isinstance(abstract_section, list):
            # Initialize an empty list to hold the strings
            abstract_parts = []
            for item in abstract_section:
                if isinstance(item, str):
                    abstract_parts.append(item)
                elif isinstance(item, dict):
                    # Assuming each dict has a 'text' key; adjust if necessary
                    abstract_parts.append(item.get('text', ''))
            # Concatenate list elements separated by a newline
            abstract = '\n'.join(abstract_parts)
        # Handle the case where abstract_section is a dict (if needed)
        elif isinstance(abstract_section, dict):
            abstract = abstract_section.get('text', '')

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

        pubmed_history = record.get('PubmedData', {}).get('History', [])

        # Extract the publication date
        article_year = None
        for item in pubmed_history:
            if item.attributes.get('PubStatus') == 'pubmed':
                article_year = item.get('Year')
                break

        venue_year = record.get('MedlineCitation', {}).get('Article', {}).get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('Year')
        if venue_year is None:
            medline_date = record.get('MedlineCitation', {}).get('Article', {}).get('Journal', {}).get('JournalIssue', {}).get('PubDate', {}).get('MedlineDate')
            if medline_date is not None:
                venue_year = medline_date.split(' ')[0]
        venue_title = record.get('MedlineCitation', {}).get('Article', {}).get('Journal', {}).get('Title')

        country = record['MedlineCitation']['MedlineJournalInfo']['Country']

        key_words = flatten_list(record.get('MedlineCitation', {}).get('KeywordList', []))

        filtered_records.append({
            'title': title,
            'abstract': abstract,
            'authors': authors,
            'pm_year': article_year,
            'country': country,
            'venue_year': venue_year,
            'venue_title': venue_title,
            'key_words': key_words,
            'url': f"https://pubmed.ncbi.nlm.nih.gov/{article_url}/"
        })

    return filtered_records


def transform_dict_to_x_y(data):
    """
    Transform the given dictionary into a format suitable for a chart.

    Args:
        data (dict): A dictionary containing the data to be transformed.

    Returns:
        list: A list of dictionaries containing the transformed data.
    """
    chart_data = []

    for key, value in data.items():
        chart_data.append({
            "x": key,
            "y": value
        })

    return chart_data


def count_entries_by_gender(data: list):
    """
    Count the number of entries in the given data by sex property.

    Args:
        data (list): A list of entries to be counted.

    Returns:
        dict: A dictionary containing the counts of entries by sex.
    """
    counts = {}

    for entry in data:
        gender = entry.get('gender')
        if gender is not None:
            if isinstance(gender, str):
                gender = [gender]
            elif not isinstance(gender, list):
                continue

            for sex_group, sex_data in sex_groups.items():
                if any(sex_value for sex_value in gender if sex_value.lower() == sex_group.lower()):
                    counts[sex_data['fda_key']] = counts.get(sex_data['fda_key'], 0) + 1

    return dict(counts.items())


def count_entries_by_age(data: list):
    """
    Count the number of entries in the given data by age property.

    Args:
        data (list): A list of entries to be counted.

    Returns:
        dict: A dictionary containing the counts of entries by age.
    """
    counts = {}

    for entry in data:
        age = entry.get('age')
        if age is not None:
            if isinstance(age, int):
                age = [age]
            elif not isinstance(age, list):
                continue

            for _, age_data in age_groups.items():
                if any(age_value for age_value in age if age_data['min'] <= age_value <= age_data['max']):
                    counts[age_data['fda_key']] = counts.get(age_data['fda_key'], 0) + 1

    return dict(counts.items())


def count_entries_by_year(data):
    """
    Count the number of entries in the given data by publication date.

    *IMPORTANT*: If the property is 'age', the algorithm counts all age groups that the entry falls into.
    For example: entry with age: [60, 120] will be put both in 'adult' and 'elderly' groups.
    This is done to avoid losing data when the age is a range since the data is super sparse.

    Args:
        data (list): A list of entries to be counted.

    Returns:
        dict: A dictionary containing the counts of entries by publication date. Sorted by year.
    """
    counts = {}

    for entry in data:
        pub_date = entry.get("pub_date")
        if pub_date:
            year = pub_date.split("-")[0]
            counts[year] = counts.get(year, 0) + 1

    return dict(sorted(counts.items(), key=lambda x: x[0]))


def count_entries_by_terms(data, search_type, limit=None):
    """
    Count the number of entries in the given data by the specified property and search criteria.

    Args:
        data (list): A list of entries to be counted.
        property (str): The property to be counted.
        search_by (str): The search criteria. Can be 'generic_name', 'brand_name', or 'side_effect'.

    Returns:
        dict: A dictionary containing the counts of entries by the specified property and search criteria. Sorted by count.
    """
    counts = {}

    for entry in data:
        if search_type == 'side_effect':
            if 'drug' in entry:
                drugs = entry['drug']
                for drug in drugs:
                    counts[drug] = counts.get(drug, 0) + 1

    if limit:
        return dict(sorted(counts.items(), key=lambda x: x[1], reverse=True)[:limit])
    else:
        return dict(sorted(counts.items(), key=lambda x: x[1], reverse=True))


def count_total_entries(data):
    """
    Count the total number of entries in the given data.

    Args:
        data (list): A list of entries to be counted.

    Returns:
        int: The total number of entries.
    """
    return len(data)


def flatten_list(data):
    """
    Flatten the given list of lists into a single list.

    Args:
        data (list): A list of lists to be flattened.

    Returns:
        list: A single list containing all the elements from the input list of lists.
    """
    return [item for sublist in data for item in sublist]