import { BackendDataType, URLParams, DemographicDataType } from '../types';
import { fetchBatchData, mapParamArrayToLabels, getNewStyleTerms } from '../utils/utils';
import { useQuery } from 'react-query';
import { backendUrl } from '../constants';

const generatePmUniversalUrls = (paramsArray: URLParams[]) => {
  const urls = paramsArray.map((params) => {
    const termsWithTypes = getNewStyleTerms(params.terms, params.searchBy);

    const url =
      `${backendUrl}/drug/get_pm_terms?` +
      `terms=${encodeURIComponent(JSON.stringify(termsWithTypes))}&` +
      `search_mode=${params.searchMode}&` +
      `sex=${params.sex}&` +
      `age=${encodeURIComponent(JSON.stringify(params.age))}&`;

    return url;
  });
  return urls;
};

const usePmDemographicData = (paramsArray: URLParams[]) => {
  const queryKey = ['pmDemographicData', JSON.stringify(paramsArray)];

  const urls = generatePmUniversalUrls(paramsArray);

  const {
    data: paramDataArray,
    error,
    isLoading,
  } = useQuery(queryKey, () => fetchBatchData(urls) as Promise<BackendDataType[]>, {
    staleTime: 3600000,
    retry: false,
    select: (data) => {
      return data.map((item, index) => {
        return {
          params: mapParamArrayToLabels(paramsArray[index] as URLParams),
          data: item.data,
          total: item.total,
        };
      }) as DemographicDataType[];
    },
    keepPreviousData: true,
  });

  const isError = !!error;

  return { paramDataArray, isError, isLoading };
};

export default usePmDemographicData;
