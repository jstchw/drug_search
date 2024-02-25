import { usePmTimeSeries } from './usePmTimeSeries';
import { useFdaTimeSeries } from './useFdaTimeSeries';

export const useDynamicTimeSeries = (isPubmed: boolean, noFilterRequest: boolean) => {
    return isPubmed ? usePmTimeSeries() : useFdaTimeSeries(noFilterRequest);
};