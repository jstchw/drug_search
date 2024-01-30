class Config:
    DEBUG = False
    FLASK_APP = 'app'
    FLASK_RUN_HOST = '0.0.0.0'
    FLASK_RUN_PORT = 16000


class DevelopmentConfig(Config):
    DEBUG = True