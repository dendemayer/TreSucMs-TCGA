summary_DF.loc[:, ('dead', case, slice(None), slice(None), slice(None))]

summary_DF.rename(columns={('dead', case, slice(None), slice(None), slice(None)):('alive', case, slice(None), slice(None), slice(None))})
import pandas as pd
pd.DataFrame().insert

Index(['UUID', 'case_id', 'gender', 'vital_status', 'drugnames',
       'survivaltime', 'years_to_last_follow_up', 'year_of_birth',
       'year_of_death', 'age_at_diagnosis', 'PROJECT', 'cut_off'],
      dtype='object')
(
['id', 'cases.0.demographic.gender', 'cases.0.demographic.vital_status',
    'pharmaceutical_therapy_drug_name', 'survivaltime', 'last_contact_days_to',
    'cases.0.demographic.year_of_birth', 'cases.0.demographic.year_of_death',
    'cases.0.diagnoses.0.age_at_diagnosis',
 'PROJECT',
 'cut_off']
Index(['survivaltime', 'cutoff', 'cases.0.demographic.gender',
       'cases.0.demographic.vital_status', 'cases.0.demographic.year_of_birth',
       'cases.0.demographic.year_of_death',
       'cases.0.diagnoses.0.age_at_diagnosis',
       'cases.0.diagnoses.0.days_to_last_follow_up', 'id',
       'bcr_patient_barcode', 'bcr_drug_barcode', 'bcr_drug_uuid',
       'form_completion_date', 'pharmaceutical_therapy_drug_name',
       'clinical_trial_drug_classification', 'pharmaceutical_therapy_type',
       'pharmaceutical_tx_started_days_to',
       'pharmaceutical_tx_ongoing_indicator',
       'pharmaceutical_tx_ended_days_to', 'treatment_best_response',
       'days_to_stem_cell_transplantation', 'pharm_regimen',
       'pharm_regimen_other', 'pharma_adjuvant_cycles_count',
       'pharma_type_other', 'pharmaceutical_tx_dose_units',
       'pharmaceutical_tx_total_dose_units', 'prescribed_dose',
       'regimen_number', 'route_of_administration',
       'stem_cell_transplantation', 'stem_cell_transplantation_type',
       'therapy_regimen', 'therapy_regimen_other', 'total_dose',
       'tx_on_clinical_trial', 'death_days_to', 'last_contact_days_to'],
      dtype='object')
pd.DataFrame().set_index()

pd.DataFrame().insert

cutoff_add_DF = cutoff_add_DF.loc[:, ['id', 'cases.0.demographic.gender', 'cases.0.demographic.vital_status', 'pharmaceutical_therapy_drug_name', 'survivaltime', 'last_contact_days_to', 'cases.0.demographic.year_of_birth', 'cases.0.demographic.year_of_death', 'cases.0.diagnoses.0.age_at_diagnosis', 'PROJECT', 'cut_off']]
