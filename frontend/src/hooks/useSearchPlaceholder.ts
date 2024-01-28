import { useEffect, useState } from "react";
import { placeholders } from "../constants";
import { SearchOptionsType } from "../types";

const useSearchPlaceholder = (
  duration: number,
  searchBy: SearchOptionsType,
) => {
  const [currentPlaceholder, setCurrentPlaceholder] = useState<
    string | undefined
  >(undefined);
  const [currentIndex, setCurrentIndex] = useState<number>(0); // used to be null

  useEffect(() => {
    const setRandomPlaceholder = () => {
      let randomIndex: number;
      do {
        randomIndex = Math.floor(
          Math.random() * (placeholders[searchBy.value]?.length || 1),
        );
      } while (randomIndex === currentIndex);

      setCurrentPlaceholder(placeholders[searchBy.value]?.[randomIndex]);
      setCurrentIndex(randomIndex);
    };

    setRandomPlaceholder(); // set a random placeholder immediately

    const interval = setInterval(setRandomPlaceholder, duration); // then update it every `duration` milliseconds

    return () => clearInterval(interval);
  }, [duration, searchBy]); // eslint-disable-line react-hooks/exhaustive-deps

  return currentPlaceholder;
};

export default useSearchPlaceholder;
