import { useUrlParams } from './useUrlParams';
import { generateFdaPath, processTermData, fetchData } from '../utils/utils';
import { useQuery } from 'react-query';
import { FDARawData, ResultItem } from 'src/types';

export const useFdaTermData = (noFilterRequest = false) => {
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

  const reportUrl = generateFdaPath(params);

  const {
    data: reportData,
    error: reportError,
    isLoading: reportLoading,
  } = useQuery(['reportUrl', reportUrl], () => fetchData(reportUrl) as Promise<FDARawData>, {
    staleTime: 3600000,
    retry: false,
    select: (data) => processTermData(data.results as ResultItem[]),
  });

  const isError = !!reportError;
  const isLoading = reportLoading;

  return { reportData, isError, isLoading };
};
