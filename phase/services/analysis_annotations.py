"""System-level analysis labels, deliberately separate from scientific states."""
from datetime import datetime, timezone
import json
import uuid

from phase.services.project_store import _write_json


def validate_scheme(payload, state_ids):
    name = str(payload.get('name') or '').strip()
    kind = payload.get('kind', 'categorical')
    categories = payload.get('categories') or []
    labels = payload.get('labels') or {}
    if not name or kind not in ('categorical', 'ordinal'):
        raise ValueError('A name and categorical/ordinal scheme kind are required.')
    if not isinstance(categories, list) or not categories or any(not isinstance(c, str) or not c.strip() for c in categories) or len(set(categories)) != len(categories):
        raise ValueError('Supply unique, nonempty categories in explicit user-defined order.')
    if not isinstance(labels, dict) or not set(labels) <= set(state_ids) or any(v not in categories for v in labels.values()):
        raise ValueError('Labels must map existing state IDs to declared categories.')
    return {'name': name, 'kind': kind, 'categories': list(categories), 'labels': dict(labels)}


def scheme_directory(store, project, system):
    store.get_system(project, system)
    return store._system_dir(project, system) / 'analysis_annotations' / 'label_schemes'


def get_scheme(store, project, system, scheme_id):
    identifier = str(uuid.UUID(str(scheme_id)))
    path = scheme_directory(store, project, system) / f'{identifier}.json'
    return json.loads(path.read_text())


def list_schemes(store, project, system):
    result = []
    for path in scheme_directory(store, project, system).glob('*.json'):
        result.append(json.loads(path.read_text()))
    return sorted(result, key=lambda s: (s['name'].casefold(), s['scheme_id']))


def save_scheme(store, project, system, payload, scheme_id=None):
    normalized = validate_scheme(payload, store.get_system(project, system).states)
    previous = get_scheme(store, project, system, scheme_id) if scheme_id else None
    identifier = previous['scheme_id'] if previous else str(uuid.uuid4())
    entry = {**normalized, 'scheme_id': identifier, 'revision': previous['revision'] + 1 if previous else 1,
             'updated_at': datetime.now(timezone.utc).isoformat()}
    _write_json(scheme_directory(store, project, system) / f'{identifier}.json', entry)
    return entry


def delete_scheme(store, project, system, scheme_id):
    entry = get_scheme(store, project, system, scheme_id)
    (scheme_directory(store, project, system) / f"{entry['scheme_id']}.json").unlink()
