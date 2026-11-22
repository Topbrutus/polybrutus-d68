from pathlib import Path
import json
R=Path(__file__).resolve().parents[1]
c=json.loads((R/'results/certificate_d68.json').read_text('utf8'))
m=json.loads((R/'publication/zenodo_metadata_blueprint.json').read_text('utf8'))
assert c['status']=='PASS' and c['dimension']==m['scientific_dimension']==68
assert c['vertices']==m['norm_one_vertices']==73
assert c['pairs']==m['exact_labelled_squared_pairs']==2628
assert c['affine_rank']==c['centered_gram_rank']==68
assert c['eigenvalue_73_over_2_multiplicity']==m['model_eigenvalue_73_over_2_multiplicity']==62
assert c['eigenvalue_36_multiplicity']==m['model_eigenvalue_36_multiplicity']==4
assert c['parent_tag_commit']=='dad596e38305066e4cb1207bb8fc947779eb445e'
assert c['parent_doi']==m['scientific_parent_doi']=='10.5281/zenodo.24432937'
assert m['related_identifiers']==[{'identifier':c['parent_doi'],'relation':'isDerivedFrom','resource_type':'publication-preprint','scheme':'doi'}]
assert len(m['keywords'])==9 and len(set(m['keywords']))==9
assert 'affine Euclidean dimension sixty-eight' in m['keywords'] and 'sixty-twofold ideal eigenvalue' in m['keywords']
assert not set(m['keywords']).intersection({'affine Euclidean dimension sixty-seven','sixty-onefold ideal eigenvalue','affine Euclidean dimension sixty-nine','sixty-threefold ideal eigenvalue'})
assert len(m['references'])==2 and c['parent_doi'] in m['references'][0] and '10.5281/zenodo.24407952' in m['references'][1]
assert all(x in (R/'publication/POLYBRUTUS_68D_NOTE.md').read_text('utf8') for x in ('## References',c['parent_doi'],'10.5281/zenodo.24407952'))
print('D68 METADATA 9 KEYWORDS / 2 REFERENCES / DIRECT D67 DOI PASS')
