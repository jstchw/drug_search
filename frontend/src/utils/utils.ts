import {
  ChartDataPoint,
  ResultItem,
  SearchTypesMap,
  TimeEventData,
  URLParams,
  SearchOptionsType,
} from "../types";
import { baseFdaUrl, searchModes, searchSex, searchTypes } from "../constants";
import React from "react";

export const processTermData = (data: ResultItem[]): ChartDataPoint[] => {
  return data.map((item) => ({
    x: item.term,
    y: item.count,
  }));
};

export const processYearData = (data: TimeEventData[]): ChartDataPoint[] => {
  const yearlyData = data.reduce(
    (acc: { [year: string]: number }, entry: TimeEventData) => {
      const year = entry.time.substring(0, 4); // Accessing the year from the time string (e.g. 2019-01-01)
      acc[year] = (acc[year] || 0) + entry.count; // If the year exists in the accumulator, add the count to it, otherwise set it to 0 and add the count to it (this is to avoid undefined errors)
      return acc;
    },
    {},
  );

  const sortedEntries: [string, number][] = Object.entries(yearlyData).sort(
    (a, b) => a[0].localeCompare(b[0]),
  );

  return sortedEntries.map(
    ([time, count]): ChartDataPoint => ({
      x: time,
      y: count,
    }),
  );
};

const formatUpperBoundDate = (date: Date) => {
  // Subtracting 1 from the year to get the previous year
  const year = date.getFullYear() - 1;

  // Constants for month and day are set to 12 and 31 respectively (last day of the year)
  const month = 12;
  const day = 31;
  return `${year}${month < 10 ? `0${month}` : month}${day < 10 ? `0${day}` : day}`;
};

export const generatePath = (
  params: URLParams,
  countType?: string,
  limitProp?: number,
) => {
  const fromDate = `20040101`;
  const limit = limitProp ? limitProp : 50;
  const toDate = formatUpperBoundDate(new Date());
  const whatToCount: SearchTypesMap = {
    generic_name: "patient.reaction.reactionmeddrapt.exact",
    brand_name: "patient.reaction.reactionmeddrapt.exact",
    receiveDate: "receivedate",
    side_effect: "patient.drug.activesubstance.activesubstancename.exact",
  };

  // Constants for min and max age
  // They are used when the age range is not fully specified
  const defaultMinAge = 0;
  const defaultMaxAge = 120;

  const searchParts = [`(receivedate:[${fromDate}+TO+${toDate}])`];

  // Term formation ----------------
  if (params.terms) {
    const encodedTerms = params.terms.flat().join("+AND+");
    searchParts.push(
      `(${mapParamToValue(params.searchBy, searchTypes)}:${encodedTerms})`,
    );
  }

  // Sex formation ----------------
  if (params.sex) {
    const sex = searchSex.find((option) => option.param === params.sex);
    if (sex) {
      searchParts.push(sex.value);
    }
  }

  // Age formation ----------------
  if (params.age) {
    if (!params.age.min) {
      params.age.min = defaultMinAge.toString();
    }

    if (!params.age.max) {
      params.age.max = defaultMaxAge.toString();
    }

    searchParts.push(
      `patient.patientonsetage:[${params.age.min}+TO+${params.age.max}]`,
    );
  }

  // Country formation ----------------
  if (params.country) {
    searchParts.push(`occurcountry:"${params.country}"`);
  }

  // countType is empty when the link is generated for Term Carousel (Feb 10, 2024)
  // countType is not empty when the link is generated for Time Series Chart (Feb 10, 2024)
  if (countType) {
    return `${baseFdaUrl}?search=${searchParts.join("+AND+")}&count=${whatToCount[countType]}`;
  } else {
    return `${baseFdaUrl}?search=${searchParts.join("+AND+")}&count=${whatToCount[params.searchBy]}&limit=${limit}`;
  }
};

export const mapParamToLabel = (
  param: string,
  options: SearchOptionsType[],
): string => {
  const option = options.find((option) => option.param === param);
  return option ? option.label : "Not Specified";
};

export const mapParamArrayToLabels = (
  params: URLParams,
): Record<string, string> => {
  const paramLabels: Record<string, string> = {};

  const assignMappedLabel = (
    param: string | null,
    optionsArray: SearchOptionsType[],
    paramName: string,
  ) => {
    const option = optionsArray.find((option) => option.param === param);
    if (option) {
      paramLabels[paramName] = option.label;
    }
  };

  if (params.searchBy) {
    assignMappedLabel(params.searchBy, searchTypes, "Type");
  }

  if (params.searchMode) {
    assignMappedLabel(params.searchMode, searchModes, "Mode");
  }

  if (params.sex) {
    assignMappedLabel(params.sex, searchSex, "Sex");
  }

  if (params.age) {
    if (params.age.max && params.age.min) {
      paramLabels["Age"] = `${params.age.min} - ${params.age.max}`;
    } else if (params.age.max && !params.age.min) {
      paramLabels["Age"] = `Up to ${params.age.max}`;
    } else if (params.age.min && !params.age.max) {
      paramLabels["Age"] = `From ${params.age.min}`;
    }
  }

  if (params.country) {
    paramLabels["Country"] = params.country;
  }

  return paramLabels;
};

export const mapParamToValue = (
  param: string,
  options: SearchOptionsType[],
): string => {
  const option = options.find((option) => option.param === param);
  return option ? option.value : "";
};

export const highlightWords = (
  text: string,
  words: string[],
): React.ReactNode => {
  const regex = new RegExp(`(${words.join("|")})`, "gi");
  const parts = text.split(regex);

  return parts.map((part, index) => {
    const isMatch = words.some((word) => new RegExp(word, "gi").test(part));
    return isMatch ? React.createElement("mark", { key: index }, part) : part;
  });
};

export const capitalizeFirstLetter = (string: string): string => {
  return string.charAt(0).toUpperCase() + string.slice(1).toLowerCase();
}