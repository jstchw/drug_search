import React, { useCallback } from "react";
import { SearchHistoryContextType, URLParams } from "../types";
import Cookies from "js-cookie";

type SearchHistoryProviderProps = {
  children: React.ReactNode;
};

export const SearchHistoryContext = React.createContext<
  SearchHistoryContextType | undefined
>(undefined);

export const SearchHistoryProvider: React.FC<SearchHistoryProviderProps> = ({
  children,
}) => {
  const initialSearchHistory = () => {
    const searchHistoryString = Cookies.get("searchHistory");
    if (searchHistoryString) {
      try {
        return JSON.parse(searchHistoryString);
      } catch (e) {
        console.error(e);
      }
    }
    return [];
  };

  const [searchHistory, setSearchHistory] =
    React.useState<URLParams[]>(initialSearchHistory);

  const updateSearchHistory = useCallback((newParams: URLParams) => {
    setSearchHistory((prevSearchHistory) => {
      const newSearchHistory = prevSearchHistory.filter(
        (item) => !compareDiffOrderStrings(item, newParams.terms),
      );

      const updatedHistory = [newParams, ...newSearchHistory].slice(0, 5); // Keep only the latest 5 entries
      Cookies.set("searchHistory", JSON.stringify(updatedHistory), {
        sameSite: "strict",
      });
      return updatedHistory;
    });
  }, []);

  const clearSearchHistory = () => {
    setSearchHistory([]);
    Cookies.remove("searchHistory");
  };

  return (
    <SearchHistoryContext.Provider
      value={{ searchHistory, updateSearchHistory, clearSearchHistory }}
    >
      {children}
    </SearchHistoryContext.Provider>
  );
};

const compareDiffOrderStrings = (item1: URLParams, item2: string[]) => {
  const sortedStr1 = item1.terms.slice().sort().join();
  const sortedStr2 = item2.slice().sort().join();
  return sortedStr1 === sortedStr2;
};
