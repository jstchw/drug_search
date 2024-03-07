import { useUrlParams } from './useUrlParams';
import { fetchData, getNewStyleTerms } from '../utils/utils';
import { useQuery } from 'react-query';
import { BackendDataType, ChartDataPoint } from 'src/types';
import { backendUrl } from '../constants';

export const usePmTermData = () => {
  let { params: { terms, searchBy, searchMode } } = useUrlParams();

  const termsWithTypes = getNewStyleTerms(terms, searchBy);

  const url =
    `${backendUrl}/drug/get_pm_terms?` +
    `terms=${encodeURIComponent(JSON.stringify(termsWithTypes))}&` +
    `search_mode=${searchMode}&`;

  const {
    data: reportData,
    error: reportError,
    isLoading: reportLoading,
  } = useQuery(['usePmTermData', url], () => fetchData(url) as Promise<BackendDataType>, {
    staleTime: 3600000,
    retry: false,
    select: (data) => data.data as ChartDataPoint[],
  });

  const isError = !!reportError;
  const isLoading = reportLoading;

  return { reportData, isError, isLoading };
};
