import { useUrlParams } from './useUrlParams';
import { fetchData, getNewStyleTerms } from '../utils/utils';
import { useQuery } from 'react-query';
import { backendUrl } from '../constants';
import { BackendDataType } from '../types';

export const usePmTimeSeries = () => {
  const {
    params: { terms, searchBy, searchMode },
  } = useUrlParams();

  const termsWithTypes = getNewStyleTerms(terms, searchBy);

  const url =
    `${backendUrl}/drug/get_pm_timedata?` +
    `terms=${encodeURIComponent(JSON.stringify(termsWithTypes))}&` +
    `search_mode=${searchMode}&`

  const {
    data,
    error: timeSeriesError,
    isLoading: timeSeriesLoading,
  } = useQuery<BackendDataType>(['pmTimeSeriesUrl', url], () => fetchData(url) as Promise<BackendDataType>, {
    staleTime: 3600000,
    retry: false,
  });

  const timeSeriesData = data?.data;
  const timeSeriesCount = data?.total;

  const isError = !!timeSeriesError;

  return { timeSeriesData, timeSeriesCount, isError, isLoading: timeSeriesLoading };
};
