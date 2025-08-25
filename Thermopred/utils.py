import os
import urllib.request

HF_BASE_URL = "https://huggingface.co/diullio/thermopred/resolve/main/"
MODELS_CACHE = os.path.join(os.path.expanduser("~"), ".thermopred", "models")

os.makedirs(MODELS_CACHE, exist_ok=True)

MODELS_URLS = {
    "MLP_gibbs.sav": HF_BASE_URL + "MLP_gibbs.sav",
    "MLP_gibbs_enthalpie.sav": HF_BASE_URL + "MLP_gibbs_enthalpie.sav",
    "RF_gibbs.sav": HF_BASE_URL + "RF_gibbs.sav",
    "RF_gibbs_Enthalpie.sav": HF_BASE_URL + "RF_gibbs_Enthalpie.sav",
    "XGB_gibbs.sav": HF_BASE_URL + "XGB_gibbs.sav",
    "XGB_gibbs_Enthalpie.sav": HF_BASE_URL + "XGB_gibbs_Enthalpie.sav",
}


def get_model(filename: str) -> str:
    """Baixa um modelo do HuggingFace se não existir em cache local"""
    local_path = os.path.join(MODELS_CACHE, filename)

    if not os.path.exists(local_path):
        print(f"[Thermopred] Baixando {filename} do HuggingFace...")
        url = MODELS_URLS[filename]
        urllib.request.urlretrieve(url, local_path)
        print(f"[Thermopred] {filename} salvo em {local_path}")

    return local_path
