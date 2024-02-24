import { URLParams } from "../types";
import { ResultItem, DemographicData } from "../types";
import {
  generatePath,
  mapParamArrayToLabels,
  processTermData,
} from "../utils/utils";
import { useQuery } from "react-query";

type UseTermDataBatchReturnType = {
  paramDataArray: DemographicData[] | undefined;
  isError: boolean;
  isLoading: boolean;
};

const fetchDemographicData = async (
  paramsArray: URLParams[],
): Promise<DemographicData[]> => {
  const fetchPromises = paramsArray.map(async (params) => {
    const url = generatePath(params, undefined, 10);
    
    const response = await fetch(url);
    if (!response.ok) throw new Error(`HTTP error: ${response.status}`);
    const result_1 = await response.json();
    return {
      params: mapParamArrayToLabels(params),
      data: processTermData(result_1.results as ResultItem[]),
    };
  });

  return Promise.all(fetchPromises);
};

export const useDemographicData = (
  paramsArray: URLParams[],
): UseTermDataBatchReturnType => {
  const queryKey = ["termDataBatch", JSON.stringify(paramsArray)];

  const {
    data: paramDataArray,
    error,
    isLoading,
  } = useQuery(queryKey, () => fetchDemographicData(paramsArray), {
    staleTime: 3600000,
    retry: false,
  });

  const isError = !!error;

  return { paramDataArray, isError, isLoading };
};
