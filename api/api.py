from flask import Flask, request, jsonify
from flask_cors import CORS

app = Flask(__name__)


@app.route("/search")
def search():
