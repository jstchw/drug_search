import { useUrlParams } from './useUrlParams';
import { useState } from 'react';
import { backendUrl } from '../constants';
import { useEffect } from 'react';
import useArticleStore from '../stores/articleStore';

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
  const { params } = useUrlParams();
  const [relevantArticles, setRelevantArticles] = useState<RelevantArticle[]>([]);
  const [relevantArticlesError, setRelevantArticlesError] = useState<unknown | boolean>(false);

  const terms = useArticleStore((state) => state.articleTerms);

  useEffect(() => {
    const fetchRelevantArticles = async () => {
      const url =
        `${backendUrl}/drug/get_articles?` +
        `terms=${params.terms}&` +
        `search_mode=${params.searchMode}&` +
        `search_type=${params.searchBy}&` +
        `sex=${params.sex}&` +
        `age=${encodeURIComponent(JSON.stringify(params.age))}&` +
        `country=${params.country}&` +
        `term_array=${encodeURIComponent(JSON.stringify(terms))}`;

      try {
        const response = await fetch(url);

        if (!response.ok) {
          setRelevantArticlesError(true);
        }

        const data = await response.json();

        setRelevantArticles(data);
      } catch (error) {
        setRelevantArticlesError(true);
        if (error instanceof Error) {
          throw new Error(error.message);
        } else {
          throw new Error('An unknown error occurred');
        }
      }
    };
    void fetchRelevantArticles();
  }, [params]);

  return { relevantArticles, relevantArticlesError };
};
