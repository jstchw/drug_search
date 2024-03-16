from sqlite3 import Date
from flask import Blueprint, jsonify, request, send_file
import io
from rdkit import Chem
from rdkit.Chem import Draw
from app.services.data_manager import DataManager
from utils.pm_utils import count_entries_by_age, count_entries_by_gender, count_entries_by_year, format_json_drug, search_json, get_pubmed_metadata, transform_dict_to_x_y, get_simple_demographic_breakdown, get_advanced_demographic_breakdown
import json
from utils.fda_utils import search_fda
import asyncio
from utils.general_utils import generate_param_list

drug_api = Blueprint('drug_api', __name__)

# DataManager initialization (getting singleton instance)
data_manager = DataManager.get_instance()
substances = data_manager.substances
products = data_manager.products
side_effects = data_manager.side_effects
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
    match search_by:
        case 'generic_name':
            return jsonify([row['Common name'] for row in substances])
        case 'brand_name':
            return jsonify([row['DrugName'] for row in products])
        case 'side_effect':
            return jsonify([row['term'] for row in side_effects])
        case _:
            return jsonify({"error": "Invalid search type"}), 400


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
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms with term value and type.
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
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms with term value and type.
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
    total_entries = len(results)

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
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms with term value and type.

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
    Retrieves the sex distribution of PubMed articles based on search parameters.

    Parameters:
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms with term value and type.

    Returns:
        A JSON response containing the total count and data of articles by sex distribution.
        If no data is found, returns a JSON response with an error message and a 404 status code.
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
        # Handle error case here
        # Handle error case here
        return jsonify({"error": "No data found"}), 404
    
    return jsonify({
        "total": sum(count_by_sex.values()),
        "data": transform_dict_to_x_y(count_by_sex)
    }), 200


@drug_api.route('/get_pm_terms', methods=['GET'])
def get_pm_terms():
    """
    Retrieves statistics about the term and the number of reports for the provided term.

    Parameters:
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms with term value and type.
    - sex (str): The sex of the patient, can be 'male' or 'female'.
    - age (int): The age of the patient.
    - country (str): The country for the search.
    - view (str): The view for the search, can be 'advanced' or 'simple'.
    - group_type (str): The group type for the search. Can be 'age' or 'sex'
    - group_name (str): The group name for the search.

    Returns:
    - JSON: A JSON response containing the total count of reports for the provided term and count by term.
    """
    params = {
        'search_mode': json.loads(request.args.get('search_mode')) if request.args.get('search_mode') else None,
        'terms': json.loads(request.args.get('terms')) if request.args.get('terms') else None,
        'sex': json.loads(request.args.get('sex')) if request.args.get('sex') else None,
        'age': json.loads(request.args.get('age')) if request.args.get('age') else None,
        'country': json.loads(request.args.get('country')) if request.args.get('country') not in [None, 'null', 'None', ''] else None,
        'view': json.loads(request.args.get('view')) if request.args.get('view') else None,
        'group_type': json.loads(request.args.get('group_type')) if request.args.get('group_type') else None,
        'return_limit': json.loads(request.args.get('return_limit')) if request.args.get('return_limit') else 10,
    }
    
    if params['view'] == 'simple':
        results = search_json(params, data=pubmed_data, limit=10000)
        demographic_breakdown = get_simple_demographic_breakdown(results, limit = params['return_limit'])
    elif params['view'] == 'advanced':
        param_list = generate_param_list(params)
        results = []
        for param in param_list:
            search_results = search_json(param, data=pubmed_data, limit=10000)
            results.extend(zip(search_results, [param] * len(search_results)))
    
        demographic_breakdown = get_advanced_demographic_breakdown(results, params['group_type'], limit = 50)

    return jsonify({
        "categories": demographic_breakdown['categories'],
        "series": demographic_breakdown['series'],
        "total_count": demographic_breakdown['total_count']
    }), 200


@drug_api.route('/get_fda_terms', methods=['GET'])
def get_fda_terms():
    """
    Retrieves statistics about the term and the number of reports for the provided term.

    Parameters:
    - search_mode (str): The search mode, can be 'relaxed' or 'strict'.
    - terms (list): A list of search terms with term value and type.
    - sex (str): The sex of the patient, can be 'male' or 'female'.
    - age (int): The age of the patient.
    - country (str): The country for the search.
    - view (str): The view for the search, can be 'advanced' or 'simple'.
    - group_type (str): The group type for the search. Can be 'age' or 'sex'
    """

    params = {
        'search_mode': json.loads(request.args.get('search_mode')) if request.args.get('search_mode') else None,
        'terms': json.loads(request.args.get('terms')) if request.args.get('terms') else None,
        'sex': json.loads(request.args.get('sex')) if request.args.get('sex') else None,
        'age': json.loads(request.args.get('age')) if request.args.get('age') else None,
        'country': json.loads(request.args.get('country')) if request.args.get('country') not in [None, 'null', 'None', ''] else None,
        'view': json.loads(request.args.get('view')) if request.args.get('view') else None,
        'group_type': json.loads(request.args.get('group_type')) if request.args.get('group_type') else None,
        'return_limit': json.loads(request.args.get('return_limit')) if request.args.get('return_limit') else 10,
    }

    print(params, flush=True)

    try:
        results = asyncio.run(search_fda(params, limit = params['return_limit']))
    except Exception as e:
        return jsonify({"error": str(e)}), 500

    return jsonify(results), 200
