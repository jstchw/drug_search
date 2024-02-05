import { useUrlParams } from "./useUrlParams";
import { useState } from "react";
import { backendUrl } from "../constants";
import { useEffect } from "react";

type RelevantArticle = {
  title: string;
  abstract: string;
  date: string;
  authors: string[];
  country: string;
  url: string;
};

export const useRelevantArticles = () => {
  const { params  } = useUrlParams();
  const [ relevantArticles, setRelevantArticles ] = useState<RelevantArticle[]>([]);
  const [ relevantArticlesError, setRelevantArticlesError ] = useState<unknown | boolean>(false);

  useEffect(() => {
    const fetchRelevantArticles = async () => {
      const url = `${backendUrl}/drug/get_articles?` +
      `term=${params.terms}&` +
      `search_type=${params.searchBy}&` +
      `sex=${params.sex}&` +
      `age=${params.age}&` +
      `country=${params.country}`;

    try {
      const response = await fetch(url);

      if (!response.ok) {
        setRelevantArticlesError(true)
      }

      const data = await response.json();

      setRelevantArticles(data);
    } catch (error) {
      setRelevantArticlesError(true);
      if (error instanceof Error) {
        throw new Error(error.message);
      } else {
        throw new Error("An unknown error occurred");
      }
    }
  }
  void fetchRelevantArticles();
  }, [params.age, params.country, params.searchBy, params.sex, params.terms]);


  return { relevantArticles, relevantArticlesError };
}