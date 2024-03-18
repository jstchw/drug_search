from flask import Flask
import logging
from flask_cors import CORS
from config import Config
from api.drug_routes import drug_api
from services.data_manager import DataManager


def create_app():
    app = Flask(__name__)
    app.config.from_object(Config)

    # Disable request logging
    log = logging.getLogger('werkzeug')
    log.setLevel(logging.ERROR)

    CORS(app, resources={r"/*": {"origins": "*"}})

    # DataManager initialization
    DataManager.get_instance()

    # Register the main API blueprint
    app.register_blueprint(drug_api, url_prefix='/api/drug')

    return app
