import { useEffect, useState } from 'react';
import { useUrlParams } from './useUrlParams';
import { optionalURLParams } from '../constants';

export const useAreParamsFiltered = () => {
  const { params } = useUrlParams();
  const [areParamsFiltered, setAreParamsFiltered] = useState<boolean>(false);

  useEffect(() => {
    const isValueSet = (value: any): boolean => {
      if (typeof value === 'object' && value !== null) {
        return Object.values(value).some(isValueSet);
      }
      return value !== null && value !== undefined;
    };

    const checkIfParamsAreFiltered = (): boolean => {
      return optionalURLParams.some((param) => {
        const value = params[param];
        if (typeof value === 'object') {
          return isValueSet(value);
        }
        return value !== null && value !== undefined;
      });
    };

    setAreParamsFiltered(checkIfParamsAreFiltered());
  }, [params]);

  return areParamsFiltered;
};
