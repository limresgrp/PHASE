import os

os.environ.setdefault("PHASE_DATA_ROOT", "/tmp/phase-reference-tests")
import pytest
pytest.importorskip("fastapi")
pytest.importorskip("httpx")
pytest.importorskip("multipart")
from fastapi import FastAPI
from fastapi.testclient import TestClient
from backend.api.v1.routes import reference, states, descriptors
from tests.test_reference_alignment import pdb, reference_system


def test_reference_upload_preview_create_and_provenance(reference_system, tmp_path, monkeypatch):
    store, _ = reference_system
    monkeypatch.setattr(reference, "project_store", store)
    monkeypatch.setattr(states, "project_store", store)
    monkeypatch.setattr(descriptors, "project_store", store)
    app = FastAPI()
    app.include_router(reference.router)
    app.include_router(states.router)
    app.include_router(descriptors.router)
    base = "/projects/p/systems/s"
    target = pdb(tmp_path / "target.pdb", "ACDEFGHIKLMNPQRSTVWY", start=108)
    with TestClient(app) as http:
        preview = http.post(base + "/reference/preview", files={"pdb": ("target.pdb", target.read_bytes())})
        assert preview.status_code == 200, preview.text
        assert preview.json()["identity"] == 1
        result = http.post(base + "/states", data={"name": "Target", "build_descriptors_after_upload": "true"}, files={"pdb": ("target.pdb", target.read_bytes())})
        assert result.status_code == 200, result.text
        created = next(iter(result.json()["states"].values()))
        assert created["alignment_status"] == "aligned"
        assert created["alignment"]["source"]["original_filename"] == "target.pdb"
        assert created["residue_keys"] == store.get_system("p", "s").reference_residue_keys
        descriptor_response = http.get(base + f"/states/{created['state_id']}/descriptors")
        assert descriptor_response.status_code == 200, descriptor_response.text
        # Target residues start at 108, but names/keys follow reference 8..27.
        assert descriptor_response.json()["residue_labels"]["res_9"] == "res_9_CYS"
        assert http.get(base + f"/states/{created['state_id']}/alignment").json()["alignment"]["quality"] == "excellent"
        assert len(http.get(base + "/alignments").json()["states"]) == 1
        assert http.get(base + "/reference").status_code == 200
        cannot_change = http.post(base + "/reference", files={"pdb": ("target.pdb", target.read_bytes())})
        assert cannot_change.status_code == 400


def test_ambiguous_api_state_does_not_register(reference_system, tmp_path, monkeypatch):
    store, ref = reference_system
    monkeypatch.setattr(states, "project_store", store)
    app = FastAPI()
    app.include_router(states.router)
    target = pdb(tmp_path / "ambiguous.pdb", "ACDEFGHIKLMNPQRSTVWY", chain="A")
    pdb(target, "ACDEFGHIKLMNPQRSTVWY", chain="B", append=True)
    with TestClient(app) as http:
        response = http.post("/projects/p/systems/s/states", data={"name": "Ambiguous"}, files={"pdb": ("target.pdb", target.read_bytes())})
        assert response.status_code == 400
        assert "Ambiguous" in response.text
        assert not store.get_system("p", "s").states


def test_empty_system_reference_endpoint_accepts_exactly_one_pdb(reference_system, monkeypatch):
    store, ref = reference_system
    store.create_system("p", name="Web", system_id="web")
    monkeypatch.setattr(reference, "project_store", store)
    monkeypatch.setattr(states, "project_store", store)
    app = FastAPI()
    app.include_router(reference.router)
    app.include_router(states.router)
    with TestClient(app) as http:
        url = "/projects/p/systems/web"
        files = [("pdb", ("one.pdb", ref.read_bytes())), ("pdb", ("two.pdb", ref.read_bytes()))]
        assert http.post(url + "/reference", files=files).status_code == 400
        response = http.post(url + "/reference", files={"pdb": ("reference.pdb", ref.read_bytes())}, data={"selection": "resid 9 to 20"})
        assert response.status_code == 200, response.text
        assert not response.json()["states"]
        first = http.post(url + "/states", data={"name": "Reference", "use_reference": "true", "build_descriptors_after_upload": "true"})
        assert first.status_code == 200, first.text
        state = next(iter(first.json()["states"].values()))
        assert state["alignment_status"] == "reference"
        assert state["n_frames"] == 1
