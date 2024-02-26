import React from 'react';
import { Button, Form, InputGroup } from 'react-bootstrap';
import { FilterLeft as FilterLeftIcon, Search as SearchIcon } from 'react-bootstrap-icons';
import './SearchBox.css';
import OptionModal from '../OptionModal/OptionModal';
import { useSuggestions } from '../../hooks/useSuggestions';
import Fuse, { FuseResult, Expression } from 'fuse.js';
import CreatableSelect from 'react-select/creatable';
import { useSearch } from '../../hooks/useSearch';
import useSearchPlaceholder from '../../hooks/useSearchPlaceholder';
import { ThemeContext } from '../../contexts/ThemeContext';
import useSearchStore from '../../stores/searchStore';

const SearchBox = () => {
  const { theme } = React.useContext(ThemeContext);

  const [searchOptions, setSearchOptions] = useSearchStore((state) => [state.searchOptions, state.setSearchOptions]);

  const [inputValue, setInputValue] = useSearchStore((state) => [state.searchInput, state.setSearchInput]);

  const inputRef = React.useRef<string[]>([]);

  const [showOptionModal, setShowOptionModal] = React.useState(false);

  const { executeSearch } = useSearch('SearchOptions');

  // Timeout for the placeholder animation
  const searchPlaceholder = useSearchPlaceholder(3000, searchOptions.searchBy);

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
      if (searchOptions.searchBy.index === 0 || searchOptions.searchBy.index === 1) {
        if (fuse) {
          const suggestion = fuse.search(e, { limit: 5 });
          if (suggestion.length > 0) {
            setSuggestions(suggestion);
          }
        }
      } else {
        /* If the search type is side effects, it will update the input value every time the user types
                       This is needed because we don't have suggestions for side effects yet (04/01/2024)
                    */
        setInputValue([]);
      }
    } else {
      setSuggestions([]);
    }
  };

  return (
    <div className={'d-flex justify-content-center text-center'}>
      <Form onSubmit={handleSearch}>
        <Form.Group controlId="searchTerm">
          <InputGroup>
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
              className={'creatable-select'}
              isMulti
              // display the passed input value if it exists
              value={inputValue.map((item) => ({ value: item, label: item }))}
              menuIsOpen={
                suggestions.length > 0 ||
                (searchOptions.searchBy.index === 2 && inputRef.current.some((item) => item.length > 0))
              }
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
