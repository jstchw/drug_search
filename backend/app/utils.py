import csv


def search_csv(column_name, search_term, file_path, return_column='DrugBank ID'):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file)

        if return_column not in reader.fieldnames or column_name not in reader.fieldnames:
            return f'Column {column_name} or {return_column} not found in the CSV file', 404

        for row in reader:
            if row[column_name].lower() == search_term.lower():
                return row[return_column]
        return 'Search term not found.', 404


def load_data(file_path, delim):
    with open(file_path, 'r') as file:
        reader = csv.DictReader(file, delimiter=delim)
        return list(reader)


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