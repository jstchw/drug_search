from sqlite3 import Date
from flask import Blueprint, jsonify, request, send_file
import io
from rdkit import Chem
from rdkit.Chem import Draw
from app.services.data_manager import DataManager
from app.utils import count_entries_by_age, count_entries_by_gender, count_entries_by_terms, count_entries_by_year, format_json_drug, search_json, get_pubmed_metadata, transform_dict_to_x_y, count_total_entries
import json
from constants import count_properties_from_search_type

drug_api = Blueprint('drug_api', __name__)

# DataManager initialization (getting singleton instance)
data_manager = DataManager.get_instance()
substances = data_manager.substances
products = data_manager.products
drug_json = data_manager.drug_json
pubmed_data = data_manager.pubmed_data


@drug_api.route('/get_molecule', methods=['GET'])
def get_molecule():
    """
    Endpoint to retrieve a molecule image based on the drug name provided in the request.
    Takes no parameters.
    Returns a molecule image as a PNG file and a status code 200 if the molecule is found, 
    otherwise returns a string "No molecule found" and a status code 404.
    """
    drug_name = request.args.get('drug_name')

    # Search the json and match the drug name to the drugbank id
    drug_id = next((row['DrugBank ID'] for row in substances if row['Common name'].lower() == drug_name.lower()), None)

    suppl = Chem.SDMolSupplier('data/molecules.sdf')

    for mol in suppl:
        if mol is not None and mol.GetProp('DRUGBANK_ID') == drug_id:
            img = Draw.MolToImage(mol)
            img_bytes = io.BytesIO()
            img.save(img_bytes, format='PNG')
            img_bytes = img_bytes.getvalue()
            return send_file(io.BytesIO(img_bytes), mimetype='image/png'), 200
    return "No molecule found", 404


@drug_api.route('/get_suggestions', methods=['GET'])
def get_suggestions():
    """
    Endpoint to retrieve real-time suggestions based on the search type provided in the request.
    Takes no parameters.
    Returns a list of suggestions as a JSON object and a status code 200 if the search type is valid, 
    otherwise returns a string "Invalid search type" and a status code 400.
    """
    search_by = request.args.get('searchBy')
    try:
        if search_by == 'patient.drug.openfda.generic_name':
            return jsonify([row['Common name'] for row in substances])
        elif search_by == 'patient.drug.openfda.brand_name':
            return jsonify([row['DrugName'] for row in products])
        else:
            return jsonify({"error": "Invalid search type"}), 400
    except Exception as e:
        return jsonify({"error": e}), 500


@drug_api.route('/get_info', methods=['GET'])
def get_info():
    """
    This function retrieves information from the drugbank json based on the provided parameters such as drug name and search type. 
    It then processes the data according to the search type and returns the filtered drug data in JSON format along with the HTTP status code 200.
    """
    # Get the params from the request
    drug_name = request.args.get('drug_name')
    search_type = request.args.get('search_type')

    # Access the drugbank json
    drug_data = drug_json['drugbank']['drug']

    if search_type == 'generic_name':
        # Search by generic name
        drug_data = [format_json_drug(drug) for drug in drug_data if drug['name'].lower() == drug_name.lower()]
        return jsonify(drug_data), 200
    elif search_type == 'brand_name':
        drug_data_filtered = []

        for drug in drug_data:
            products = drug.get('products', {}).get('product')

            # Check if 'products' is a single dictionary or a list of dictionaries
            if isinstance(products, dict):
                products = [products]  # Convert to a list for uniform handling

            # Now 'products' is always a list
            if isinstance(products, list):
                if any(product.get('name', '').lower() == drug_name.lower() for product in products):
                    drug_data_filtered.append(format_json_drug(drug, drug_name.capitalize()))

        return jsonify(drug_data_filtered), 200


@drug_api.route('/get_articles', methods=['GET'])
def get_articles():
    """
    Retrieves articles based on the provided search parameters and returns a JSON response.

    Parameters:
    - search_type (str): The type of search, can be 'generic_name', 'brand_name', or 'side_effect'.
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms.
    - sex (str): The gender, can be 'male' or 'female'.
    - age (int): The age of the person.
    - country (str): The country for the search.

    Returns:
    - JSON: A JSON response containing the articles based on the search parameters.
    """
    params = {
        'search_mode': request.args.get('search_mode') if request.args.get('search_mode') in ['relaxed', 'strict'] else None,
        'sex': request.args.get('sex') if request.args.get('sex') in ['male', 'female'] else None,
        'age': json.loads(request.args.get('age')) if request.args.get('age') else None,
        'country': request.args.get('country') if request.args.get('country') not in [None, 'null', 'None', ''] else None,
        'terms': json.loads(request.args.get('terms')) if request.args.get('terms') else None
    }

    if not params['terms'] or not params['search_mode']:
        return jsonify({"error": "Missing required parameters"}), 400

    results = search_json(params,
                          data=pubmed_data,
                          limit=10)

    # Form a list of pubmed ids and use it to generate the pubmed url
    if results:
        pubmed_ids = [result['pubmed_id'] for result in results]

        article_metadata = get_pubmed_metadata(pubmed_ids)

        return jsonify(article_metadata), 200
    else:
        return jsonify({"error": "No articles found"}), 404


@drug_api.route('/get_pm_timedata', methods=['GET'])
def get_pm_timedata():
    """
    Retrieves number of publications for the provided term and returns a ready-to-use
    dictionary in format {total: 'total count of publications', [x: 'year', y: 'number of publications'...]}.

    Parameters:
    - search_type (str): The type of search, can be 'generic_name', 'brand_name', or 'side_effect'.
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms.
    - sex (str): The sex of the patient, can be 'male' and 'female'.
    - age (int): The age of the patient.
    - country (str): The country for the search.

    Returns:
    - JSON: A JSON response containing the total count of publications for the provided term and count by year.
    """
    params = {
        'search_mode': request.args.get('search_mode') if request.args.get('search_mode') in ['relaxed', 'strict'] else None,
        'terms': json.loads(request.args.get('terms')) if request.args.get('terms') else None,
        'sex': request.args.get('sex') if request.args.get('sex') in ['male', 'female'] else None,
        'age': json.loads(request.args.get('age')) if request.args.get('age') else None,
        'country': request.args.get('country') if request.args.get('country') not in [None, 'null', 'None', ''] else None
    }

    results = search_json(params, data=pubmed_data, limit=10000)
    total_entries = count_total_entries(results)

    if results:
        count_by_year = count_entries_by_year(results)
    else:
        count_by_year = {"error": "No data found"}

    if not results or not count_by_year:
        return jsonify({"error": "No data found"}), 404

    # Remove years outside the range of start_year and end_year
    start_year = 2004
    end_year = Date.today().year - 1
    count_by_year = {year: count for year, count in count_by_year.items() if start_year <= int(year) <= end_year}

    return jsonify({
        "total": total_entries,
        "data": transform_dict_to_x_y(count_by_year)
    }), 200


@drug_api.route('/get_pm_age_distribution', methods=['GET'])
def get_pm_age_distribution():
    """
    Retrieves number of publications for the provided term and returns a ready-to-use
    dictionary in format {total: count, data: [x: 'age_group', y: 'number of publications'...]}.

    Parameters:
    - search_type (str): The type of search, can be 'generic_name', 'brand_name', or 'side_effect'.
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms.

    Returns:
    - JSON: A JSON response containing the total count of publications for the provided term and count by age.
    """

    params = {
        'search_mode': request.args.get('search_mode') if request.args.get('search_mode') in ['relaxed', 'strict'] else None,
        'terms': json.loads(request.args.get('terms')) if request.args.get('terms') else None,
    }

    results = search_json(params, data=pubmed_data, limit=1000)

    if results:
        count_by_age = count_entries_by_age(results)
    else:
        count_by_age = {"error": "No data found"}

    if not results or not count_by_age:
        return jsonify({"error": "No data found"}), 404

    return jsonify({
        "total": sum(count_by_age.values()),
        "data": transform_dict_to_x_y(count_by_age)
    }), 200


@drug_api.route('/get_pm_sex_distribution', methods=['GET'])
def get_pm_sex_distribution():
    """
    Retrieves number of publications for the provided params and returns a ready-to-use
    dictionary in format {total: count, data: [x: 'age_group', y: 'number of publications'...]}.

    Parameters:
    - search_type (str): The type of search, can be 'generic_name', 'brand_name', or 'side_effect'.
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms.

    Returns:
    - JSON: A JSON response containing the total count of publications for the provided term and count by age.
    """
    
    params = {
        'search_mode': request.args.get('search_mode') if request.args.get('search_mode') in ['relaxed', 'strict'] else None,
        'terms': json.loads(request.args.get('terms')) if request.args.get('terms') else None,
    }

    results = search_json(params, data=pubmed_data, limit=1000)

    if results:
        count_by_sex = count_entries_by_gender(results)
    else:
        count_by_sex = {"error": "No data found"}

    if not results or not count_by_sex:
        return jsonify({"error": "No data found"}), 404
    
    return jsonify({
        "total": sum(count_by_sex.values()),
        "data": transform_dict_to_x_y(count_by_sex)
    }), 200


@drug_api.route('/get_pm_terms', methods=['GET'])
def get_pm_terms():
    """
    Retrieves statistics about the term and an amount of reports for the provided term.

    Parameters:

    Returns:
    - JSON: Ready-to-use dictionary in format {total: count, data: [x: 'term', y: 'number of reports'...]}.
    """
    params = {
        'search_mode': request.args.get('search_mode') if request.args.get('search_mode') in ['relaxed', 'strict'] else None,
        'terms': json.loads(request.args.get('terms')) if request.args.get('terms') else None,
        'sex': request.args.get('sex') if request.args.get('sex') in ['male', 'female'] else None,
        'age': json.loads(request.args.get('age')) if request.args.get('age') else None,
        'country': request.args.get('country') if request.args.get('country') not in [None, 'null', 'None', ''] else None
    }

    results = search_json(params, data=pubmed_data, limit=1000)
    
    # Current solution to count by term works only if all terms are of the same type
    # It can be changed to accommodate different way of searching but this will be done when
    # an AI model will be able to determine the type of each term automatically.
    if results:
        count_by_term = count_entries_by_terms(results, params['terms'][0]['type'], limit=50)
    else:
        count_by_term = {"error": "No data found"}

    if not results or not count_by_term:
        return jsonify({"error": "No data found"}), 404


    return jsonify({
        "total": sum(count_by_term.values()),
        "data": transform_dict_to_x_y(count_by_term)
    }), 200