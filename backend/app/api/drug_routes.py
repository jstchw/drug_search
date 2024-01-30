from flask import Blueprint, jsonify, request, send_file
import io
from rdkit import Chem
from rdkit.Chem import Draw
from app.utils import search_csv
from app.services.data_manager import DataManager
from app.utils import format_json_drug

drug_api = Blueprint('drug_api', __name__)

# DataManager initialization (getting singleton instance)
data_manager = DataManager.get_instance()
substances = data_manager.substances
products = data_manager.products
drug_json = data_manager.drug_json


@drug_api.route('/get_molecule', methods=['GET'])
def get_molecule():
    drug_name = request.args.get('drug_name')
    drug_id = search_csv('Common name', drug_name, 'data/suggestions/suggestion_db.csv')
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
