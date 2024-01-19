export type SearchOptions = {
    searchBy: SearchOptionsType;
    sex: SearchOptionsType;
    age: AgeOptions;
    country: SearchOptionsType;
}

export type SearchOptionsType = {
    value: string;
    index: number;
    label: string;
    type: string;
    enabled?: boolean;
    param?: string;
}

export type AgeOptions = {
    enabled: boolean;
    min: SearchOptionsType;
    max: SearchOptionsType;
}

export type DrugProperties = {
    name: string;
    classification?: string;
    groups?: string[];
    iupac?: string;
    formula?: string;
    brands?: string[];
    half_life?: string;
    indication?: string;
    product?: string;
    molecule_url?: string;
    [key: string]: string | string[] | undefined;
}

export type URLParams = {
    terms: string[];
    searchBy: string;
    sex?: string | null;
    age?: {
        min: string | null;
        max: string | null;
    };
    country?: string | null;
}

export type SearchHistoryContextType = {
    searchHistory: URLParams[];
    updateSearchHistory: (params: URLParams) => void;
    clearSearchHistory: () => void;
}

export type SearchTypesMap = {
    generic_name: string;
    brand_name: string;
    receiveDate: string;
    side_effect: string;
    [key: string]: string;
}

export type PlaceholderType = {
    [key: string]: string[];
}

export type Results = {
    [key: string]: number;
}

export type ResultItem = {
    term: string;
    count: number;
}

export type FDARawData = {
    meta:
    {
        disclaimer: string;
        terms: string;
        license: string;
        last_updated: string;
    }
    results: ResultItem[] | TimeEventData[];
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

export type DemographicGroups = {
    name: string;
    sex: string;
    age: [number, number];
    def: string;
}