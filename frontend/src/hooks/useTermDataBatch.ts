import { FDARawData, URLParams } from "../types";
import { ChartDataPoint, ResultItem } from "../types";
import {
  generatePath,
  mapParamArrayToLabels,
  processTermData,
} from "../utils/utils";
import { useQuery } from "react-query";

type BatchData = {
  params: Record<string, string>;
  data: ChartDataPoint[];
};

type UseTermDataBatchReturnType = {
  paramDataArray: BatchData[] | undefined;
  error: boolean;
  isLoading: boolean;
};

const fetchBatchData = async (
  paramsArray: URLParams[],
): Promise<BatchData[]> => {
  const fetchPromises = paramsArray.map((params) => {
    const url = generatePath(params, undefined, 10);
    return fetch(url)
      .then((response) => {
        if (!response.ok) throw new Error(`HTTP error: ${response.status}`);
        return response.json();
      })
      .then((result: FDARawData) => ({
        params: mapParamArrayToLabels(params),
        data: processTermData(result.results as ResultItem[]),
      }));
  });

  return Promise.all(fetchPromises);
};

export const useTermDataBatch = (
  paramsArray: URLParams[],
): UseTermDataBatchReturnType => {
  const queryKey = ["termDataBatch", JSON.stringify(paramsArray)];

  const {
    data: paramDataArray,
    error,
    isLoading,
  } = useQuery(queryKey, () => fetchBatchData(paramsArray), {
    staleTime: 3600000,
    retry: false,
  });

  const isError = !!error;

  return { paramDataArray, error: isError, isLoading };
};
