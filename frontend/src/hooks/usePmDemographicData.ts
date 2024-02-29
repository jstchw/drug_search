import { URLParams } from "../types";
import { fetchBatchData } from "../utils/utils";
import { useQuery } from "react-query";
import { backendUrl } from "../constants";

const generatePmUniversalUrls = (paramsArray: URLParams[]) => {
    const urls = paramsArray.map((params) => {
        const url = `${backendUrl}/drug/get_pm_terms?` +
            `terms=${params.terms}&` +
            `search_mode=${params.searchMode}&` +
            `search_type=${params.searchBy}&` +
            `sex=${params.sex}&` +
            `age_min=${params.age?.min}&` +
            `age_max=${params.age?.max}&`;

        return url;
    });
    return urls;
}

const usePmDemographicData = (paramsArray: URLParams[]) => {
    const queryKey = ['pmDemographicData', JSON.stringify(paramsArray)];

    const urls = generatePmUniversalUrls(paramsArray);

    const {
        data: paramDataArray,
        error,
        isLoading,
    } = useQuery(queryKey, () => fetchBatchData(urls), {
        staleTime: 3600000,
        retry: false,
    });

    const isError = !!error;

    return { paramDataArray, isError, isLoading };
    
}

export default usePmDemographicData;