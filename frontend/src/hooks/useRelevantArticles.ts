import { useUrlParams } from './useUrlParams';
import { backendUrl } from '../constants';
import { useQuery } from 'react-query';
import useArticleStore from '../stores/articleStore';
import { fetchData } from '../utils/utils';

type RelevantArticle = {
  title: string;
  abstract: string;
  pm_year: string;
  authors: string[];
  country: string;
  venue_title: string;
  venue_year: string;
  key_words: string[];
  url: string;
};

export const useRelevantArticles = () => {
  const terms = useArticleStore((state) => state.articleTerms);
  const { params } = useUrlParams();

  const url =
    `${backendUrl}/drug/get_articles?` +
    `search_mode=${params.searchMode}&` +
    `sex=${params.sex}&` +
    `age=${encodeURIComponent(JSON.stringify(params.age))}&` +
    `country=${params.country}&` +
    `terms=${encodeURIComponent(JSON.stringify(terms))}`;

  const {
    data = [],
    error,
    isLoading,
    isFetching,
  } = useQuery<RelevantArticle[]>(
    ['relevantArticlesUrl', url],
    () => {
      if (terms.length > 0) {
        return fetchData(url) as Promise<RelevantArticle[]>;
      } else {
        return Promise.resolve([]);
      }
    },
    {
      staleTime: 3600000,
      retry: false,
      keepPreviousData: true,
      refetchOnWindowFocus: false,
    }
  );

  const loading = isLoading || isFetching;

  return { relevantArticles: data, relevantArticlesError: error, isLoading: loading };

  return { relevantArticles: data, relevantArticlesError: error, isLoading };
};
