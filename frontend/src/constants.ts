// Search types that can be selected in the popover
import countries from './assets/countries.json';
import {
  PlaceholderType,
  SearchOptionsType,
  AgeOptions,
  AgeGroups,
  DrugGroupConfig,
  PercentageIntensityColors,
  SearchTypesMap,
} from './types';
import { Prescription, MagnifyingGlass, PawPrint, Trash, Carrot, Flask, X, SmileyNervous } from '@phosphor-icons/react';

export const searchTypes: SearchOptionsType[] = [
  {
    value: 'patient.drug.openfda.generic_name',
    index: 0,
    label: 'Generic Name',
    type: 'searchBy',
    param: 'generic_name',
  },
  {
    value: 'patient.drug.openfda.brand_name',
    index: 1,
    label: 'Brand Name',
    type: 'searchBy',
    param: 'brand_name',
  },
  {
    value: 'patient.reaction.reactionmeddrapt',
    index: 2,
    label: 'Side Effect',
    type: 'searchBy',
    param: 'side_effect',
  },
];

export const searchModes: SearchOptionsType[] = [
  {
    value: 'relaxed',
    index: 0,
    label: 'Relaxed',
    type: 'searchMode',
    param: 'relaxed',
  },
  {
    value: 'strict',
    index: 1,
    label: 'Strict',
    type: 'searchMode',
    param: 'strict',
  },
];
export const searchSex: SearchOptionsType[] = [
  {
    value: 'patient.patientsex:1',
    index: 0,
    label: 'Male',
    type: 'sex',
    param: 'male',
  },
  {
    value: 'patient.patientsex:2',
    index: 1,
    label: 'Female',
    type: 'sex',
    param: 'female',
  },
];

export const searchAgeRange: AgeOptions = {
  enabled: false,
  ageGroupsEnabled: true,
  min: {
    value: '0',
    index: 0,
    label: 'Minimum Age',
    type: 'age',
    param: 'min_age',
  },
  max: {
    value: '120',
    index: 1,
    label: 'Maximum Age',
    type: 'age',
    param: 'max_age',
  },
};

export const searchAgeGroups: AgeGroups = {
  '0 - 2': {
    min: 0,
    max: 2,
  },
  '3 - 11': {
    min: 3,
    max: 11,
  },
  '12 - 17': {
    min: 12,
    max: 17,
  },
  '18 - 64': {
    min: 18,
    max: 64,
  },
  '65 - 85': {
    min: 65,
    max: 85,
  },
  '85+': {
    min: 85,
    max: 120,
  },
};

export const ageGroupsFDA: Record<string, string> = {
  1: 'Neonate',
  2: 'Infant',
  3: 'Child',
  4: 'Adolescent',
  5: 'Adult',
  6: 'Elderly',
};

export const sexGroupsFDA: Record<string, string> = {
  0: 'Unknown',
  1: 'Male',
  2: 'Female',
};

export const whatToCount: SearchTypesMap = {
  generic_name: 'patient.reaction.reactionmeddrapt.exact',
  brand_name: 'patient.reaction.reactionmeddrapt.exact',
  receiveDate: 'receivedate',
  side_effect: 'patient.drug.activesubstance.activesubstancename.exact',
  age_group: 'patient.patientagegroup',
  patient_sex: 'patient.patientsex',
};

export const searchCountry: SearchOptionsType[] = Object.entries(countries).map(([value, label], index) => ({
  value,
  index,
  label,
  type: 'country',
  param: value,
}));

export const defaultSearchOptions = {
  searchBy: { ...(searchTypes[0] as SearchOptionsType), enabled: true },
  searchMode: { ...searchModes[0], enabled: true },
  sex: searchSex[0],
  age: searchAgeRange,
  country: searchCountry.find(({ value }) => value === 'US') as SearchOptionsType,
};

export const optionalURLParams = ['sex', 'age', 'country'];

// redo the placeholder data structure
export const placeholders: PlaceholderType = {
  [searchTypes[0]?.value || 'patient.drug.openfda.generic_name']: [
    'Acetaminophen',
    'Alprazolam',
    'Oxycodone',
    'Modafinil',
    'Ibuprofen',
    'Diazepam',
    'Hydrocodone',
    'Tramadol',
    'Codeine',
    'Gabapentin',
    'Meloxicam',
    'Cyclobenzaprine',
    'Naproxen',
    'Methocarbamol',
    'Prednisone',
    'Citalopram',
    'Amitriptyline',
    'Trazodone',
    'Lisinopril',
    'Atorvastatin',
    'Metformin',
    'Amlodipine',
    'Omeprazole',
    'Metoprolol',
    'Simvastatin',
    'Losartan',
    'Azithromycin',
    'Hydrochlorothiazide',
    'Amoxicillin',
    'Albuterol',
    'Levothyroxine',
    'Furosemide',
    'Fluticasone',
    'Montelukast',
    'Fluoxetine',
    'Escitalopram',
    'Sertraline',
    'Bupropion',
    'Cyclobenzaprine',
    'Trazodone',
    'Lisinopril',
    'Atorvastatin',
    'Metformin',
    'Amlodipine',
    'Omeprazole',
    'Metoprolol',
    'Simvastatin',
    'Losartan',
    'Azithromycin',
    'Hydrochlorothiazide',
    'Amoxicillin',
    'Albuterol',
    'Levothyroxine',
    'Furosemide',
    'Fluticasone',
    'Montelukast',
    'Fluoxetine',
    'Escitalopram',
    'Sertraline',
    'Bupropion',
  ],
  [searchTypes[1]?.value || 'patient.drug.openfda.brand_name']: [
    'Xanax',
    'Percocet',
    'Adderall',
    'Valium',
    'Vicodin',
    'Ambien',
    'Klonopin',
    'Oxycontin',
    'Concerta',
    'Ritalin',
    'Zoloft',
    'Ativan',
    'Lyrica',
    'Lunesta',
    'Lexapro',
    'Prozac',
  ],
  [searchTypes[2]?.value || 'patient.reaction.reactionmeddrapt']: ['Headache', 'Nausea', 'Fever'],
};

export const drugGroupConfig: DrugGroupConfig = {
  approved: {
    label: 'Approved',
    variant: 'success',
    description: 'Has been approved in at least one jurisdiction, at some point in time.',
    IconComponent: Prescription,
  },
  investigational: {
    label: 'Investigational',
    variant: 'warning',
    description: 'Undergoing evaluation in the drug approval process in at least one jurisdiction.',
    IconComponent: MagnifyingGlass,
  },
  illicit: {
    label: 'Illicit',
    variant: 'danger',
    description: 'Deemed illegal in at least one jurisdiction, at some point in time.',
    IconComponent: X,
  },
  vet_approved: {
    label: 'Vet Approved',
    variant: 'info',
    description: 'Approved for veterinary use in at least one jurisdiction, at some point in time.',
    IconComponent: PawPrint,
  },
  withdrawn: {
    label: 'Withdrawn',
    variant: 'secondary',
    description: 'Has been withdrawn from the market in at least one jurisdiction, at some point in time.',
    IconComponent: Trash,
  },
  nutraceutical: {
    label: 'Nutraceutical',
    variant: 'light',
    description: 'Pharmaceutical-grade nutrient with potential health benefits.',
    IconComponent: Carrot,
  },
  experimental: {
    label: 'Experimental',
    variant: 'dark',
    description: 'Shown to bind specific proteins in experimental settings.',
    IconComponent: Flask,
  },
  side_effect: {
    label: 'Side Effect',
    variant: 'primary',
    description: 'Has been reported as a side effect.',
    IconComponent: SmileyNervous,
  },
};

export const backendUrl = `${window.location.protocol}//${window.location.hostname}:16000/api`;

export const analysisFrontendUrl = 'http://drug_watch_analysis_frontend:3001';

export const baseFdaUrl = 'https://api.fda.gov/drug/event.json';

export const fuseOptions = {
  keys: ['name'],
  threshold: 0.4,
  minMatchCharLength: 2,
};

export const versionInfo = {
  appName: 'DrugWatch',
  tag: 'Alpha',
  number: '0.1',
};

export const percentageIntensityColors: PercentageIntensityColors = {
  extremelyCommon: {
    percentageRange: [50, 100],
    color: '#B71C1C', // Deep Red
  },
  veryCommon: {
    percentageRange: [20, 50],
    color: '#C62828', // Strong Red
  },
  moderatelyCommon: {
    percentageRange: [10, 20],
    color: '#D32F2F', // Muted Red
  },
  common: {
    percentageRange: [5, 10],
    color: '#FFA000', // Amber
  },
  lessCommon: {
    percentageRange: [4, 5],
    color: '#FFB300', // Slightly darker Amber
  },
  common3: {
    percentageRange: [3, 4],
    color: '#FFC107', // Standard Amber
  },
  common2: {
    percentageRange: [2, 3],
    color: '#AED581', // Light Green
  },
  common1: {
    percentageRange: [1, 2],
    color: '#81C784', // Medium Light Green
  },
  uncommon: {
    percentageRange: [0.5, 1],
    color: '#4DB6AC', // Teal Green
  },
  rare: {
    percentageRange: [0.1, 0.5],
    color: '#26A69A', // Darker Teal
  },
  veryRare: {
    percentageRange: [0.01, 0.1],
    color: '#00897B', // Deep Teal
  },
  extremelyRare: {
    percentageRange: [0, 0.01],
    color: '#757575', // Grey
  },
};

export const chartColors = [
  '#59768A',
  '#035363',
  '#32B2BF',
  '#D5E0BE',
  '#CE9062',
  '#E0AB86',
  '#C7CE8A',
  '#6EB585',
  '#325951',
  '#6F9F9D',
];

export const oppositeAggregation: Record<string, string> = {
  Age: 'Sex',
  Sex: 'Age',
};
