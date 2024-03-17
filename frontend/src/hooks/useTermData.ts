import { useUrlParams } from './useUrlParams';
import { fetchData, getNewStyleTerms } from '../utils/utils';
import { useQuery } from 'react-query';
import { DemographicResponseType } from 'src/types';
import { backendUrl, searchAgeRange } from '../constants';

export const useTermData = (source: string, noFilterRequest = false) => {
  const { params } = useUrlParams();

  const termsWithTypes = getNewStyleTerms(params.terms, params.searchBy);

  const baseUrl = `${backendUrl}/drug/get_${source}_terms?`;
  const termsParam = `terms=${encodeURIComponent(JSON.stringify(termsWithTypes))}`;
  const searchModeParam = `search_mode=${encodeURIComponent(JSON.stringify(params.searchMode))}`;
  const sexParam = params.sex ? `&sex=${encodeURIComponent(JSON.stringify(params.sex))}` : '';
  const viewParam = `&view=${encodeURIComponent(JSON.stringify('simple'))}`;
  const return_limit = `&return_limit=${encodeURIComponent(JSON.stringify(50))}`;

  const ageParam = noFilterRequest
    ? `&age=${encodeURIComponent(JSON.stringify({
        min: searchAgeRange.min.value,
        max: searchAgeRange.max.value,
      }))}`
    : params.age
    ? `&age=${encodeURIComponent(JSON.stringify(params.age))}`
    : '';

  let url: string;
  if (noFilterRequest) {
    url = baseUrl + termsParam + '&' + searchModeParam + viewParam + ageParam + return_limit;
  } else {
    url = baseUrl + termsParam + '&' + searchModeParam + sexParam + ageParam + viewParam + return_limit;
  }

  const { data, isError, isLoading } = useQuery(
    ['usePmTermData', url],
    () => fetchData(url) as Promise<DemographicResponseType>,
    {
      staleTime: 3600000,
      retry: false,
    }
  );

  return { data, isError, isLoading };
};