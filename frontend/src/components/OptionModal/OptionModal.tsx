import React, { useEffect } from "react";
import {
  Button,
  FigureCaption,
  Form,
  InputGroup,
  Modal,
  Nav,
  ToggleButton,
  Tooltip,
  Overlay,
} from "react-bootstrap";
import { ThemeContext } from "../../contexts/ThemeContext";
import {
  Moon,
  Sun,
  Funnel,
  ClockCounterClockwise,
} from "@phosphor-icons/react";
import Cookies from "js-cookie";
import {
  searchCountry,
  searchSex,
  searchTypes,
  searchModes,
} from "../../constants";
import { SearchOptions, SearchOptionsType } from "../../types";
import SearchHistory from "../SearchHistory/SearchHistory";
import { Carousel } from "react-responsive-carousel";

const OptionModal = (props: {
  showOptionModal: boolean;
  setShowOptionModal: React.Dispatch<boolean>;
  searchOptions: SearchOptions;
  setSearchOptions: React.Dispatch<SearchOptions>;
}) => {
  useEffect(() => {
    const optionsString = JSON.stringify(props.searchOptions);
    Cookies.set("searchOptions", optionsString, {
      expires: 365,
      sameSite: "strict",
    });
  }, [props.searchOptions]);

  const { theme, toggleTheme } = React.useContext(ThemeContext);

  const [carouselIndex, setCarouselIndex] = React.useState(0);

  const ageTooltipTarget = React.useRef(null);
  const [incorrectAge, setIncorrectAge] = React.useState(false);

  const minAgeRef = React.useRef<HTMLInputElement>(null);
  const maxAgeRef = React.useRef<HTMLInputElement>(null);

  const incorrectAgeTooltip = (
    <Tooltip id="incorrectAgeTooltip">
      <span>Invalid age. Please enter a valid number or range.</span>
    </Tooltip>
  );

  const handleSearchTypeChange = (value: string) => {
    const newSearchType =
      (searchTypes.find(
        (searchType) => searchType.value === value,
      ) as SearchOptionsType) || searchTypes[0];

    props.setSearchOptions({
      ...props.searchOptions,
      searchBy: {
        ...newSearchType,
        enabled: props.searchOptions.searchBy.enabled ?? true,
      },
    });
  };

  const handleSearchModeChange = (value: string) => {
    const newSearchMode =
      (searchModes.find(
        (searchMode) => searchMode.value === value,
      ) as SearchOptionsType) || searchModes[0];

    props.setSearchOptions({
      ...props.searchOptions,
      searchMode: {
        ...newSearchMode,
        enabled: props.searchOptions.searchMode.enabled ?? true,
      },
    });
  };

  const handleSexChange = (value: string) => {
    const newSearchSex =
      (searchSex.find(
        (searchSex) => searchSex.value === value,
      ) as SearchOptionsType) || searchSex[0];

    props.setSearchOptions({
      ...props.searchOptions,
      sex: {
        ...newSearchSex,
        enabled: props.searchOptions.sex.enabled ?? true,
      },
    });
  };

  const handleAgeChange = (e: React.ChangeEvent<HTMLInputElement>) => {
    const isNumeric = /^\d*$/.test(e.currentTarget.value);

    if (
      minAgeRef.current &&
      maxAgeRef.current &&
      parseInt(minAgeRef.current.value) > parseInt(maxAgeRef.current.value)
    ) {
      setIncorrectAge(true);
    } else {
      setIncorrectAge(false);
    }

    if (isNumeric) {
      const index = parseInt(e.currentTarget.id);

      const updatedAge = { ...props.searchOptions.age };

      if (index === updatedAge.min.index) {
        updatedAge.min = { ...updatedAge.min, value: e.currentTarget.value };
      } else if (index === updatedAge.max.index) {
        updatedAge.max = { ...updatedAge.max, value: e.currentTarget.value };
      }

      props.setSearchOptions({
        ...props.searchOptions,
        age: updatedAge,
      });
    }
  };

  const handleCountryChange = (value: string) => {
    const newSearchCountry =
      (searchCountry.find(
        (searchCountry) => searchCountry.index === parseInt(value),
      ) as SearchOptionsType) || searchCountry[0];

    props.setSearchOptions({
      ...props.searchOptions,
      country: {
        ...newSearchCountry,
        enabled: props.searchOptions.country.enabled ?? true,
      },
    });
  };

  return (
    <>
      <Modal
        centered
        show={props.showOptionModal}
        onHide={() => !incorrectAge && props.setShowOptionModal(false)}
      >
        <Modal.Header closeButton>
          {/*Theme change button*/}
          <Button
            variant={"outline-primary"}
            onClick={toggleTheme}
            className={"d-flex align-items-center"}
          >
            {theme === "light" ? <Moon /> : <Sun />}
          </Button>
          <div className={"vr mx-2"} />
          <Modal.Title>Options</Modal.Title>
        </Modal.Header>
        <Modal.Body>
          <Carousel
            showThumbs={false}
            showIndicators={false}
            showArrows={false}
            showStatus={false}
            selectedItem={carouselIndex}
            dynamicHeight={true}
          >
            <Form>
              {/* Search type change */}
              <Form.Group className="mb-3 d-flex align-items-center">
                <div className={"d-flex align-items-center"}>
                  <ToggleButton
                    id={"search_type_button"}
                    type="checkbox"
                    variant="outline-primary"
                    checked={true}
                    value="1"
                    disabled={true}
                  >
                    Search by
                  </ToggleButton>
                </div>
                <InputGroup
                  className={"flex-grow-1 mx-3"}
                  style={{ width: "auto" }}
                >
                  <Form.Select
                    onChange={(e) =>
                      handleSearchTypeChange(e.currentTarget.value)
                    }
                    value={props.searchOptions.searchBy.value}
                    style={{ width: "auto" }}
                  >
                    {searchTypes.map((searchType, index) => (
                      <option key={index} value={searchType.value}>
                        {searchType.label}
                      </option>
                    ))}
                  </Form.Select>
                </InputGroup>
              </Form.Group>

              {/* Search mode change */}
              <Form.Group className="mb-3 d-flex align-items-center">
                <div className={"d-flex align-items-center"}>
                  <ToggleButton
                    id={"search_mode_button"}
                    type="checkbox"
                    variant="outline-primary"
                    checked={true}
                    value="1"
                    disabled={true}
                  >
                    Search mode
                  </ToggleButton>
                </div>
                <InputGroup
                  className={"flex-grow-1 mx-3"}
                  style={{ width: "auto" }}
                >
                  <Form.Select
                    onChange={(e) =>
                      handleSearchModeChange(e.currentTarget.value)
                    }
                    value={props.searchOptions.searchMode.value}
                    style={{ width: "auto" }}
                  >
                    {searchModes.map((searchMode, index) => (
                      <option key={index} value={searchMode.value}>
                        {searchMode.label}
                      </option>
                    ))}
                  </Form.Select>
                </InputGroup>
              </Form.Group>

              {/* Sex option change */}
              <Form.Group className="mb-3 d-flex align-items-center">
                <div className={"d-flex align-items-center"}>
                  <ToggleButton
                    id={"sex_change_button"}
                    type="checkbox"
                    variant="outline-primary"
                    checked={props.searchOptions.sex.enabled}
                    value="1"
                    onClick={() =>
                      props.setSearchOptions({
                        ...props.searchOptions,
                        sex: {
                          ...(props.searchOptions.sex as SearchOptionsType),
                          enabled: !props.searchOptions.sex.enabled,
                        },
                      })
                    }
                  >
                    Sex
                  </ToggleButton>
                </div>
                <InputGroup className={"mx-3 flex-grow-1"}>
                  <Form.Select
                    onChange={(e) => {
                      handleSexChange(e.currentTarget.value);
                    }}
                    value={props.searchOptions.sex.value}
                    disabled={!props.searchOptions.sex.enabled}
                  >
                    {searchSex.map((sex, index) => (
                      <option key={index} value={sex.value}>
                        {sex.label}
                      </option>
                    ))}
                  </Form.Select>
                </InputGroup>
              </Form.Group>

              {/* Age option change */}
              <Form.Group className="mb-3 d-flex align-items-center">
                <div className={"d-flex align-items-center"}>
                  <Overlay
                    placement={"left"}
                    show={incorrectAge}
                    target={ageTooltipTarget.current}
                  >
                    {incorrectAgeTooltip}
                  </Overlay>
                  <ToggleButton
                    ref={ageTooltipTarget}
                    id={"age_change_button"}
                    type="checkbox"
                    variant="outline-primary"
                    checked={props.searchOptions.age.enabled}
                    value="1"
                    onClick={() => {
                      const newSearchOptions = { ...props.searchOptions };
                      newSearchOptions.age = {
                        ...props.searchOptions.age,
                        enabled: !props.searchOptions.age.enabled,
                      };
                      props.setSearchOptions(newSearchOptions);
                    }}
                  >
                    Age
                  </ToggleButton>
                </div>

                {/* Min age input */}
                <InputGroup className={"mx-3 flex-grow-1"}>
                  <Form.Control
                    ref={minAgeRef}
                    type="text"
                    placeholder="Min"
                    value={props.searchOptions.age.min.value ?? 0}
                    id={"0"}
                    onChange={(e) => {
                      handleAgeChange(e as React.ChangeEvent<HTMLInputElement>);
                    }}
                    disabled={!props.searchOptions.age.enabled}
                    onCopy={(e) => e.preventDefault()}
                    onPaste={(e) => e.preventDefault()}
                  />

                  <InputGroup.Text>-</InputGroup.Text>

                  {/* Max age input */}
                  <Form.Control
                    ref={maxAgeRef}
                    type="text"
                    placeholder="Max"
                    value={props.searchOptions.age.max.value ?? 0}
                    id={"1"}
                    onChange={(e) => {
                      handleAgeChange(e as React.ChangeEvent<HTMLInputElement>);
                    }}
                    disabled={!props.searchOptions.age.enabled}
                    onCopy={(e) => e.preventDefault()}
                    onPaste={(e) => e.preventDefault()}
                  />
                </InputGroup>
              </Form.Group>

              {/* Country option change */}
              <Form.Group>
                <div className={"d-flex align-items-center"}>
                  <ToggleButton
                    id={"country_change_button"}
                    type="checkbox"
                    variant="outline-primary"
                    value="1"
                    onClick={() =>
                      props.setSearchOptions({
                        ...props.searchOptions,
                        country: {
                          ...(props.searchOptions.country as SearchOptionsType),
                          enabled: !props.searchOptions.country.enabled,
                        },
                      })
                    }
                    checked={props.searchOptions.country.enabled}
                  >
                    Country
                  </ToggleButton>
                  <InputGroup className={"mx-3 flex-grow-1"}>
                    <Form.Select
                      disabled={!props.searchOptions.country.enabled}
                      onChange={(e) =>
                        handleCountryChange(e.currentTarget.value)
                      }
                      value={props.searchOptions.country.index}
                    >
                      {searchCountry.map((country, index) => (
                        <option key={index} value={index}>
                          {country.label}
                        </option>
                      ))}
                    </Form.Select>
                  </InputGroup>
                </div>
              </Form.Group>
            </Form>

            <SearchHistory setShowOptionModal={props.setShowOptionModal} />
          </Carousel>

          <Nav
            variant="tabs"
            defaultActiveKey={carouselIndex}
            className={"mt-3"}
          >
            <Nav.Item>
              <Nav.Link
                className={"d-flex align-items-center"}
                eventKey="0"
                onClick={() => setCarouselIndex(0)}
              >
                <Funnel weight={"light"} />
                <div className={"vr mx-2"} />
                Filters
              </Nav.Link>
            </Nav.Item>
            <Nav.Item>
              <Nav.Link
                className={"d-flex align-items-center"}
                eventKey="1"
                onClick={() => setCarouselIndex(1)}
              >
                <ClockCounterClockwise weight={"light"} />
                <div className={"vr mx-2"} />
                History
              </Nav.Link>
            </Nav.Item>
          </Nav>

          <FigureCaption className={"d-flex justify-content-center mt-3"}>
            <a
              className={"text-decoration-none"}
              href={"https://github.com/jstchw/drug_search"}
              target={"_blank"}
              rel="noreferrer"
            >
              DrugSearch Alpha 0.1.1
            </a>
          </FigureCaption>
        </Modal.Body>
      </Modal>
    </>
  );
};

export default OptionModal;
