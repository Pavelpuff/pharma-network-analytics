import uvicorn
from fastapi import FastAPI, HTTPException
from pydantic import BaseModel

from src.data_parser import DrugBankParser

# Application initialization
app = FastAPI(
    title="Pharma Analytics API",
    description="API for retrieving pharmacological pathways and drug interactions.",
    version="1.0.0",
)

# Load data into memory at startup using the professional parser module
PARSER = DrugBankParser("data/drugbank_partial.xml")
try:
    DF_PATHWAYS = PARSER.extract_pathways()
except Exception as e:
    DF_PATHWAYS = None
    print(f"Warning: Could not load XML database. {e}")


class DrugRequest(BaseModel):
    drug_name: str


class DrugResponse(BaseModel):
    drug_name: str
    interaction_count: int


@app.post("/drug_interactions/", response_model=DrugResponse)
def get_drug_interactions(request: DrugRequest):
    """
    Returns the number of pathways with which a given drug interacts.
    """
    if DF_PATHWAYS is None:
        raise HTTPException(status_code=500, detail="Database not loaded properly.")

    count = DF_PATHWAYS[DF_PATHWAYS["Drug_Name"] == request.drug_name].shape[0]

    if count == 0:
        raise HTTPException(
            status_code=404,
            detail=f"Drug '{request.drug_name}' not found or has no pathways.",
        )

    return DrugResponse(drug_name=request.drug_name, interaction_count=count)


if __name__ == "__main__":
    uvicorn.run("src.api:app", host="0.0.0.0", port=8000, reload=True)
