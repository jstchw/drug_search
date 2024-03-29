import React, { useEffect } from 'react';
import { Button, Form, InputGroup } from 'react-bootstrap';
import { FilterLeft as FilterLeftIcon, Search as SearchIcon } from 'react-bootstrap-icons';
import OptionModal from '../OptionModal/OptionModal';
import { useSuggestions } from '../../hooks/useSuggestions';
import Fuse, { FuseResult, Expression } from 'fuse.js';
import CreatableSelect from 'react-select/creatable';
import { useSearch } from '../../hooks/useSearch';
import useSearchPlaceholder from '../../hooks/useSearchPlaceholder';
import useSearchStore from '../../stores/searchStore';
import useGeneralOptionsStore from '../../stores/generalOptionsStore';
import { CSSObjectWithLabel } from 'react-select';

const SearchBox = () => {
  const theme = useGeneralOptionsStore((state) => state.theme);

  const [searchOptions, setSearchOptions] = useSearchStore((state) => [state.searchOptions, state.setSearchOptions]);

  const [inputValue, setInputValue] = useSearchStore((state) => [state.searchInput, state.setSearchInput]);

  const inputRef = React.useRef<string[]>([]);
  const [isInputRefPresent, setIsInputRefPresent] = React.useState(false);

  const [showOptionModal, setShowOptionModal] = React.useState(false);

  const { executeSearch } = useSearch('SearchOptions');

  // Timeout for the placeholder animation
  const searchPlaceholder = useSearchPlaceholder(3000, searchOptions.searchBy, isInputRefPresent);

  const [fuse, setFuse] = React.useState<Fuse<string> | null>(null);
  // Elements for the suggestion mechanism
  const [suggestions, setSuggestions] = React.useState<FuseResult<string>[]>([]);

  useSuggestions(searchOptions.searchBy, setFuse);

  const handleSearch = async (e: { preventDefault: () => void }) => {
    if (e) e.preventDefault();
    executeSearch(inputValue, searchOptions);
  };

  const handleInputChange = (e: string | Expression) => {
    // inputRef is used for tracking the input value in real time without state delays
    // Needed to manipulate the menu (open and close) when there are no suggestions
    inputRef.current = [e.toString()];
    if (e) {
      if (fuse) {
        const suggestion = fuse.search(e, { limit: 5 });
        if (suggestion.length > 0) {
          setSuggestions(suggestion);
        }
      }
    } else {
      setSuggestions([]);
    }
  };

  // Stop the placeholder animation when the input is not empty
  useEffect(() => {
    setIsInputRefPresent(inputRef.current.some((element) => element.trim() !== ''));
  }, [inputRef.current]);

  return (
    <div className={'d-flex justify-content-center text-center'}>
      <Form onSubmit={handleSearch}>
        <Form.Group controlId="searchTerm">
          <InputGroup className="flex-nowrap">
            <>
              <Button
                variant="outline-secondary"
                id="button-filters"
                onClick={() => setShowOptionModal(true)}
                className={'d-flex align-items-center'}
              >
                <FilterLeftIcon />
              </Button>
              <OptionModal
                showOptionModal={showOptionModal}
                setShowOptionModal={setShowOptionModal}
                searchOptions={searchOptions}
                setSearchOptions={setSearchOptions}
              />
            </>

            <CreatableSelect
              isMulti
              // display the passed input value if it exists
              value={inputValue.map((item) => ({ value: item, label: item }))}
              menuIsOpen={suggestions.length > 0 || inputRef.current.some((item) => item.length > 0)}
              onInputChange={(e) => handleInputChange(e)}
              onChange={(e) => setInputValue(e.map((item) => item.value))}
              options={suggestions.map((suggestion) => ({
                value: suggestion.item,
                label: suggestion.item,
              }))}
              components={{
                DropdownIndicator: () => null,
                IndicatorSeparator: () => null,
              }}
              placeholder={searchPlaceholder}
              classNames={{
                control: () =>
                  theme === 'dark'
                    ? 'bg-transparent border-secondary rounded-0'
                    : 'bg-white border-secondary rounded-0',
                menu: () => (theme === 'dark' ? 'bg-dark text-light' : ''),
                input: () => (theme === 'dark' ? 'bg-transparent text-light' : ''),
                option: (state) => (state.isFocused && theme === 'dark' ? 'text-dark' : ''),
                multiValue: () => 'bg-secondary-subtle',
                multiValueLabel: () => (theme === 'dark' ? 'text-light' : ''),
              }}
              styles={{
                control: (provided) =>
                  ({
                    ...provided,
                    minWidth: '11rem',
                  }) as CSSObjectWithLabel,
              }}
            />

            <Button variant="outline-primary" id="button-submit" type="submit" className={'d-flex align-items-center'}>
              <SearchIcon />
            </Button>
          </InputGroup>
        </Form.Group>
      </Form>
    </div>
  );
};

export default SearchBox;
