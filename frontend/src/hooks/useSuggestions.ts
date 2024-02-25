import React, { useEffect } from 'react';
import Fuse from 'fuse.js';
import { backendUrl, fuseOptions } from '../constants';
import { SearchOptionsType } from '../types';

export const useSuggestions = (searchBy: SearchOptionsType, setFuse: React.Dispatch<Fuse<string>>) => {
  useEffect(() => {
    async function fetchSuggestions() {
      try {
        if (searchBy.param !== 'side_effect') {
          // Form a URL based on the searchBy parameter and fetch the data
          const url = `${backendUrl}/drug/get_suggestions?searchBy=${encodeURIComponent(searchBy.value)}`;
          const response = await fetch(url);
          const data: string[] = await response.json();
          setFuse(new Fuse(data, fuseOptions));
        }
      } catch (error) {
        console.error('Error fetching suggestions:', error);
      }
    }

    fetchSuggestions();
  }, [searchBy, setFuse]);
};
