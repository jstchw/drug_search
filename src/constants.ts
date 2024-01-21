// Search types that can be selected in the popover
import countries from "./assets/countries.json";
import {PlaceholderType, SearchOptionsType, AgeOptions} from "./types";

export const searchTypes: SearchOptionsType[] = [
    {
        value: 'patient.drug.openfda.generic_name',
        index: 0,
        label: 'Generic Name',
        type: 'searchBy',
        param: 'generic_name'
    },
    {
        value: 'patient.drug.openfda.brand_name',
        index: 1,
        label: 'Brand Name',
        type: 'searchBy',
        param: 'brand_name'
    },
    {
        value: 'patient.reaction.reactionmeddrapt',
        index: 2,
        label: 'Side Effect',
        type: 'searchBy',
        param: 'side_effect'
    },
]
export const searchSex: SearchOptionsType[] = [
    {
        value: 'patient.patientsex:1',
        index: 0,
        label: 'Male',
        type: 'sex',
        param: 'male'
    },
    {
        value: 'patient.patientsex:2',
        index: 1,
        label: 'Female',
        type: 'sex',
        param: 'female'
    }
]

export const searchAgeRange: AgeOptions = {
    enabled: false,
    min: {
        value: '0',
        index: 0,
        label: 'Minimum Age',
        type: 'age',
        param: 'min_age'
    },
    max: {
        value: '100',
        index: 1,
        label: 'Maximum Age',
        type: 'age',
        param: 'max_age'
    },
}
export const searchCountry: SearchOptionsType[] = Object.entries(countries).map(([value, label], index) => ({
    value,
    index,
    label,
    type: 'country',
    param: value
}))

export const defaultSearchOptions = {
    searchBy: { ...searchTypes[0] as SearchOptionsType, enabled: true },
    sex: searchSex[0],
    age: searchAgeRange,
    country: searchCountry.find(({ value }) => value === 'US') as SearchOptionsType,
}

// redo the placeholder data structure
export const placeholders: PlaceholderType = {
    [searchTypes[0]?.value || 'patient.drug.openfda.generic_name']: ['Acetaminophen', 'Alprazolam', 'Oxycodone', 'Modafinil', 'Ibuprofen', 'Diazepam',
        'Hydrocodone', 'Tramadol', 'Codeine', 'Gabapentin', 'Meloxicam', 'Cyclobenzaprine', 'Naproxen',
        'Methocarbamol', 'Prednisone', 'Citalopram', 'Amitriptyline', 'Trazodone', 'Lisinopril', 'Atorvastatin',
        'Metformin', 'Amlodipine', 'Omeprazole', 'Metoprolol', 'Simvastatin', 'Losartan', 'Azithromycin',
        'Hydrochlorothiazide', 'Amoxicillin', 'Albuterol', 'Levothyroxine', 'Furosemide', 'Fluticasone',
        'Montelukast', 'Fluoxetine', 'Escitalopram', 'Sertraline', 'Bupropion', 'Cyclobenzaprine', 'Trazodone',
        'Lisinopril', 'Atorvastatin', 'Metformin', 'Amlodipine', 'Omeprazole', 'Metoprolol', 'Simvastatin', 'Losartan',
        'Azithromycin', 'Hydrochlorothiazide', 'Amoxicillin', 'Albuterol', 'Levothyroxine', 'Furosemide',
        'Fluticasone', 'Montelukast', 'Fluoxetine', 'Escitalopram', 'Sertraline', 'Bupropion'],
    [searchTypes[1]?.value || 'patient.drug.openfda.brand_name']: ['Xanax', 'Percocet', 'Adderall', 'Valium', 'Vicodin', 'Ambien', 'Klonopin', 'Oxycontin',
        'Concerta', 'Ritalin', 'Zoloft', 'Ativan', 'Lyrica', 'Lunesta', 'Lexapro', 'Prozac'],
    [searchTypes[2]?.value || 'patient.reaction.reactionmeddrapt']: ['Headache', 'Nausea', 'Fever'],
};

export const drugGroups = {
    'approved': {
        description: 'Has been approved in at least one jurisdiction, at some point in time.',
        badgeVariant: 'success'
    },
    'investigational': {
        description: 'Undergoing evaluation in the drug approval process in at least one jurisdiction.',
        badgeVariant: 'warning'
    },
    'illicit': {
        description: 'Deemed illegal in at least one jurisdiction, at some point in time.',
        badgeVariant: 'danger'
    },
    'vet_approved': {
        description: 'Approved for veterinary use in at least one jurisdiction, at some point in time.',
        badgeVariant: 'info'
    },
    'withdrawn': {
        description: 'Has been withdrawn from the market in at least one jurisdiction, at some point in time.',
        badgeVariant: 'secondary'
    },
    'nutraceutical': {
        description: 'Pharmaceutical-grade nutrient with potential health benefits.',
        badgeVariant: 'light'
    },
    'experimental': {
        description: 'Shown to bind specific proteins in experimental settings.',
        badgeVariant: 'dark'
    }
}

export const backendUrl = process.env.NODE_ENV === 'development' ?
    `${window.location.protocol}//${window.location.hostname}:16000/api/` :
    'https://drugsearch.org/api/'

export const baseFdaUrl = 'https://api.fda.gov/drug/event.json'

export const fuseOptions = {
    keys: ['name'],
    threshold: 0.4,
    minMatchCharLength: 2,
}