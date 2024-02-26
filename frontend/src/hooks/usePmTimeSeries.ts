import { useUrlParams } from './useUrlParams';
import { fetchData } from '../utils/utils';
import { useQuery } from 'react-query';
import { backendUrl } from '../constants';
import { BackendDataType } from '../types';

export const usePmTimeSeries = () => {
  const {
    params: { terms, searchBy, searchMode },
  } = useUrlParams();

  const url =
    `${backendUrl}/drug/get_pm_timedata?` +
    `terms=${terms}&` +
    `search_mode=${searchMode}&` +
    `search_type=${searchBy}&`;

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
