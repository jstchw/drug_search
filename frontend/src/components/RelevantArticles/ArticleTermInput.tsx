import CreatableSelect from 'react-select/creatable';
import { CSSObjectWithLabel } from 'react-select';
import { useState, useRef} from 'react';
import { useSuggestions } from '../../hooks/useSuggestions';
import Fuse, { FuseResult, Expression } from 'fuse.js';
import { SearchOptionsType } from '../../types';
import { searchTypes } from '../../constants';
import { Form, Dropdown, InputGroup } from 'react-bootstrap';
import useArticleStore from '../../stores/articleStore';
import useGeneralOptionsStore from '../../stores/generalOptionsStore';
import { motion } from 'framer-motion';

const ArticleTermInput = () => {
  const theme = useGeneralOptionsStore((state) => state.theme);

  const [inputValue, setInputValue] = useState<string>('');
  const [searchType, setSearchType] = useState<SearchOptionsType>(searchTypes[0]!);
  const [fuse, setFuse] = useState<Fuse<string> | null>(null);
  const [suggestions, setSuggestions] = useState<FuseResult<string>[]>([]);
  const inputRef = useRef<string[]>([]);

  const addArticleTerm = useArticleStore((state) => state.addArticleTerm);

  const [isDropdownOpen, setIsDropdownOpen] = useState(false);
  const handleDropdownToggle = () => setIsDropdownOpen(!isDropdownOpen);

  useSuggestions(searchType, setFuse);

  const handleInputChange = (e: string | Expression) => {
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

  const handleSubmit = (value: string) => {
    addArticleTerm(value, searchType.param, true);
    setInputValue('');
  }

  return (
    <motion.div layout>
      <Form className={'mb-2'}>
        <InputGroup>
          <Dropdown onToggle={handleDropdownToggle} show={isDropdownOpen}>
            <Dropdown.Toggle 
              variant="outline-secondary" 
              id="dropdown-basic"
            >
              {searchType.label}
            </Dropdown.Toggle>


            <Dropdown.Menu as={motion.div}
              initial="closed"
              animate={isDropdownOpen ? "open" : "closed"}
              variants={{
                  open: { opacity: 1 },
                  closed: { opacity: 0 }
              }}
              transition={{ duration: 0.2, ease: "easeInOut" }}
            >
              {searchTypes.map((type, index) => (
                <Dropdown.Item
                  key={index}
                  onClick={() => setSearchType(type)}
                >
                  {type.label}
                </Dropdown.Item>
              ))}
            </Dropdown.Menu>
          </Dropdown>

          <CreatableSelect
            className={'text-center'}
            menuIsOpen={
              suggestions.length > 0 ||
              inputRef.current.some((item) => item.length > 0)
            }
            components={{
              DropdownIndicator: () => null,
              IndicatorSeparator: () => null,
            }}
            onInputChange={(e) => handleInputChange(e)}
            options={suggestions.map((suggestion) => ({
              value: suggestion.item,
              label: suggestion.item,
            }))}
            onChange={(e) => e && handleSubmit(e.value)}
            value={ inputValue ?
              { value: inputValue, label: inputValue } : null
            }
            placeholder={'Filter by terms...'}
            classNames={{
              control: () => 
              theme === 'dark'
                    ? 'bg-transparent border-secondary rounded-start-0 rounded-end-2'
                    : 'bg-white border-secondary rounded-start-0 rounded-end-2',
                    menu: () => (theme === 'dark' ? 'bg-dark text-light' : ''),
                    input: () => (theme === 'dark' ? 'bg-transparent text-light' : ''),
                    option: (state) => (state.isFocused && theme === 'dark' ? 'text-dark' : ''),
            }}
            styles={{
              control: (base) => ({
                ...base,
                minWidth: '13em',
              }) as CSSObjectWithLabel,
            }}
          />
        </InputGroup>
      </Form>
    </motion.div>
  );
}

export default ArticleTermInput;