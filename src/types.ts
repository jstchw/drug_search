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

export type SearchHistoryElement = {
    terms: string[];
    options: SearchOptions;
}