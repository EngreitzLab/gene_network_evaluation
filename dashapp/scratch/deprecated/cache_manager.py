from flask_caching import Cache
from dash import Dash

def setup_cache(app: Dash):
    # Setup the cache using Flask-Caching and attach it to the app
    cache = Cache(app.server, config={
        'CACHE_TYPE': 'simple',  # Adjust the cache type as needed
        'CACHE_DEFAULT_TIMEOUT': 300
    })
    return cache