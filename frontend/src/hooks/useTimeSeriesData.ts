import { generatePath } from '../utils/utils';
import { useUrlParams } from './useUrlParams';
import { processYearData } from '../utils/utils';
import { fetchData } from '../utils/utils';
import { useQuery } from 'react-query';
import { TimeEventData } from 'src/types';

export const useTimeSeriesData = (noFilterRequest = false) => {
  let { params } = useUrlParams();

  if (noFilterRequest) {
    params = {
      ...params,
      sex: null,
      age: {
        min: null,
        max: null,
      },
      country: null,
    };
  }

  const timeSeriesUrl = generatePath(params, 'receiveDate');

  const seriesCountUrl = generatePath(params, undefined, undefined, true);

  const {
    data: timeSeriesData,
    error: timeSeriesError,
    isLoading: timeSeriesLoading,
  } = useQuery(['timeSeriesUrl', timeSeriesUrl], () => fetchData(timeSeriesUrl), {
    staleTime: 3600000,
    retry: false,
    select: (data) => processYearData(data.results as TimeEventData[]),
  });

  const {
    data: timeSeriesCount,
    error: seriesCountError,
    isLoading: seriesCountLoading,
  } = useQuery(['seriesCountUrl', seriesCountUrl], () => fetchData(seriesCountUrl), {
    staleTime: 3600000,
    retry: false,
    select: (data) => data.meta.results && data.meta.results.total,
  });

  const isError = !!timeSeriesError || !!seriesCountError;
  const isLoading = timeSeriesLoading || seriesCountLoading;

  return { timeSeriesData, timeSeriesCount, isError, isLoading };
};
