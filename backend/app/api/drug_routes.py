from flask import Blueprint, jsonify, request, send_file
import io
from rdkit import Chem
from rdkit.Chem import Draw
from app.services.data_manager import DataManager
from app.utils import format_json_drug, search_json, get_pubmed_metadata

drug_api = Blueprint('drug_api', __name__)

# DataManager initialization (getting singleton instance)
data_manager = DataManager.get_instance()
substances = data_manager.substances
products = data_manager.products
drug_json = data_manager.drug_json


@drug_api.route('/get_molecule', methods=['GET'])
def get_molecule():
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
    params = {
        'search_type': request.args.get('search_type'),
        'terms': request.args.get('term').strip().split(','),
        'sex': request.args.get('sex'),
        'age': request.args.get('age'),
        'country': request.args.get('country')
    }

    if not params['search_type'] or not params['terms']:
        return jsonify({"error": "Missing required parameters"}), 400

    results = search_json(params,
                          json_file_path='data/pubmed_article_data/pubmed_data.json',
                          limit=10)

    # Form a list of pubmed ids and use it to generate the pubmed url
    if results:
        pubmed_ids = [result['pubmed_id'] for result in results]

        article_metadata = get_pubmed_metadata(pubmed_ids)

        return jsonify(article_metadata), 200
    else:
        return jsonify({"error": "No articles found"}), 404