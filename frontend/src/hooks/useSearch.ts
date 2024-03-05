import { useNavigate } from 'react-router-dom';
import { SearchOptions, AgeOptions, SearchOptionsType, URLParams } from '../types';
import useArticleStore from '../stores/articleStore';

export const useSearch = (paramType: 'SearchOptions' | 'URLParams') => {
  const navigate = useNavigate();
  const executeSearch = (inputValue: string[], params: SearchOptions | URLParams) => {
    if (inputValue.length > 0) {
      const prettyInputValues = inputValue.map((value) => value.charAt(0).toLowerCase() + value.slice(1));
      const queryParams = new URLSearchParams();
      queryParams.append('query', prettyInputValues.join('-'));

      // Reset the article matching term on new search
      useArticleStore.getState().resetArticleTerms();

      if (paramType === 'SearchOptions') {
        const searchOptions = params as SearchOptions;
        Object.entries(searchOptions).forEach(([key, value]) => {
          if (value !== null && value !== undefined) {
            if (key === 'age' && value.enabled) {
              const ageOptions = value as AgeOptions;
              if (ageOptions.min.value) {
                queryParams.append('min_age', ageOptions.min.value);
              }
              if (ageOptions.max.value) {
                queryParams.append('max_age', ageOptions.max.value);
              }
            } else if (key !== 'age' && value.enabled) {
              const searchParams = value as SearchOptionsType;
              queryParams.append(key, searchParams.param || searchParams.value);
            }
          }
        });
      } else if (paramType === 'URLParams') {
        const urlParams = params as URLParams;
        Object.entries(urlParams).forEach(([key, value]) => {
          if (value !== null && value !== undefined) {
            if (key === 'age') {
              if (params.age?.min) {
                queryParams.append('min_age', params.age.min as string);
              }
              if (params.age?.max) {
                queryParams.append('max_age', params.age.max as string);
              }
            } else if (key !== 'age' && key !== 'terms') {
              queryParams.append(key, typeof value === 'string' ? value : JSON.stringify(value));
            }
          }
        });
      }

      navigate(`/search?${queryParams.toString()}`);
    }
  };

  return { executeSearch };
};
