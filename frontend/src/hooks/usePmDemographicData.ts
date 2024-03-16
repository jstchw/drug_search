import { BackendDataType, URLParams, DemographicDataType, DemographicRequestType } from '../types';
import { fetchData, mapParamArrayToLabels, getNewStyleTerms } from '../utils/utils';
import { useQuery } from 'react-query';
import { backendUrl } from '../constants';

interface PmDemographicData {
  categories: string[];
  series: number[];
  total_count: number;
}

const usePmDemographicData = (requestArgs: DemographicRequestType) => {
  const queryKey = ['pmDemographicData', JSON.stringify(requestArgs)];

  const baseUrl = `${backendUrl}/drug/get_pm_terms?`;
  const termsParam = `terms=${encodeURIComponent(JSON.stringify(requestArgs.terms))}`;
  const searchModeParam = `search_mode=${encodeURIComponent(JSON.stringify(requestArgs.searchMode))}`;
  const sexParam = requestArgs.sex ? `&sex=${encodeURIComponent(JSON.stringify(requestArgs.sex))}` : '';
  const ageParam = requestArgs.age ? `&age=${encodeURIComponent(JSON.stringify(requestArgs.age))}` : '';
  const groupTypeParam = `&group_type=${encodeURIComponent(JSON.stringify(requestArgs.groupType))}`;
  const viewParam = `&view=${encodeURIComponent(JSON.stringify(requestArgs.view))}`;

  const url = baseUrl + termsParam + '&' + searchModeParam + sexParam + ageParam + groupTypeParam + viewParam;

  const {
    data,
    error,
    isLoading,
  } = useQuery(queryKey, () => fetchData(url) as Promise<PmDemographicData>, {
    staleTime: 3600000,
    retry: false,
    keepPreviousData: true,
  });

  const isError = !!error;

  return { data, isError, isLoading };
};

export default usePmDemographicData;
