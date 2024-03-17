import { DemographicRequestType, DemographicResponseType } from '../types';
import { useQuery } from 'react-query';
import { backendUrl } from '../constants';
import { fetchData } from '../utils/utils';

const useDemographicData = (requestArgs: DemographicRequestType, source: string) => {
  const baseUrl = `${backendUrl}/drug/get_${source}_terms?`;
  const termsParam = `terms=${encodeURIComponent(JSON.stringify(requestArgs.terms))}`;
  const searchModeParam = `search_mode=${encodeURIComponent(JSON.stringify(requestArgs.searchMode))}`;
  const sexParam = requestArgs.sex ? `&sex=${encodeURIComponent(JSON.stringify(requestArgs.sex))}` : '';
  const ageParam = requestArgs.age ? `&age=${encodeURIComponent(JSON.stringify(requestArgs.age))}` : '';
  const groupTypeParam = `&group_type=${encodeURIComponent(JSON.stringify(requestArgs.groupType))}`;
  const viewParam = `&view=${encodeURIComponent(JSON.stringify(requestArgs.view))}`;

  const url = baseUrl + termsParam + '&' + searchModeParam + sexParam + ageParam + groupTypeParam + viewParam;

  const queryKey = ['useDemographicData' + JSON.stringify(requestArgs), source];

  const { data, error, isLoading, isFetching } = useQuery(
    queryKey,
    () => fetchData(url) as Promise<DemographicResponseType>,
    {
      staleTime: 3600000,
      retry: false,
      keepPreviousData: true,
    }
  );

  const isError = !!error;

  const loading = isLoading || isFetching;

  return { data, isError, isLoading: loading };
};

export default useDemographicData;
