export type SearchOptions = {
    searchBy: string;
    sex: string | undefined;
    age: {
        value: string;
        index: number;
        type: string;
    } | undefined;
    country: string | undefined;
}

export type SearchOptionsType = {
    value: string;
    index: number;
    label: string;
    type: string;
}

export type SearchHistoryElement = {
    terms: string[];
    options: SearchOptions;
}

export type Results = {
    [key: string]: number;
}

export type TimeEventData = {
    time: string;
    count: number;
}

export type ChartDataPoint = {
    x: string;
    y: number;
}

export type ThemeType = 'light' | 'dark';

export type DrugInfo = {
    drug_name: string;
    drug_class: string;
    groups: string[];
    half_life: string;
    indication: string;
    iupac: string;
    molecule_url: string;
    product: string;
    brands: string[];
    ADE: null;
    formula: string;
    full_info: boolean;
}