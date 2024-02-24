import { useUrlParams } from "../../hooks/useUrlParams";
import { useNavigate } from "react-router-dom";
import SearchBox from "../../components/SearchBox/SearchBox";
import Header from "../../components/Header/Header";
import { Col, Row } from "react-bootstrap";
import useDrugInfo from "../../hooks/useDrugInfo";
import React from "react";
import DrugPropertyBox from "../../components/DrugPropertyBox/DrugPropertyBox";
import { SearchHistoryContext } from "../../contexts/SearchHistoryContext";
import ChartSection from "../../components/ChartSection/ChartSection";
import { DrugProperties } from "../../types";
import { PropagateLoader } from "react-spinners";
import { Bug } from "@phosphor-icons/react";
import { isMobile } from "react-device-detect";
import RelevantArticles from "../../components/RelevantArticles/RelevantArticles";
import InfoCard from "../../components/InfoCard/InfoCard";
import Scrollbar from "../../components/Scrollbar/Scrollbar";
import useSearchStore from "../../stores/searchStore";
import { AnimatePresence, useInView } from "framer-motion";

const Result = () => {
  const { params, paramError } = useUrlParams();

  const setSearchInput = useSearchStore(state => state.setSearchInput);

  const capitalizedTerms = React.useMemo(
    () =>
      params.terms.map((term) => term.charAt(0).toUpperCase() + term.slice(1)),
    [params.terms], // Only recalculate when params.terms changes
  );

  React.useEffect(() => {
    if (capitalizedTerms.length > 0) {
      setSearchInput(capitalizedTerms);
    }
  }, [capitalizedTerms, setSearchInput]);

  const navigate = useNavigate();

  const [loading, setLoading] = React.useState(true);

  const primaryColor = window
    .getComputedStyle(document.documentElement)
    .getPropertyValue("--primary")
    .trim();

  // Retrieve drug info from the API (DrugSearch server)
  const { drugInfo, drugInfoError } = useDrugInfo(params);

  React.useEffect(() => {
    setLoading(true);
  }, [params]);

  React.useEffect(() => {
    if (drugInfo.length > 0 || drugInfoError) {
      setLoading(false);
    }
  }, [drugInfo, drugInfoError]);

  // Grouping drugs by brand name only if the search type is brand name (specified in the conditional statement during rendering)
  const groupedByBrandName = React.useMemo(() => {
    return drugInfo.reduce<Record<string, DrugProperties[]>>((acc, drug) => {
      if (drug.product) {
        (acc[drug.product] = acc[drug.product] || []).push(drug);
      }
      return acc;
    }, {});
  }, [drugInfo]);

  // Effect to redirect to the error page if there is an error
  React.useEffect(() => {
    if (paramError) {
      navigate("/error");
    }
  }, [paramError, navigate]);

  // Effect to update the document title
  React.useEffect(() => {
    if (params.terms && capitalizedTerms.length > 0) {
      document.title = capitalizedTerms.join(" & ") + " - Drug Search";
    }
  }, [capitalizedTerms, params.terms]);

  // Updating search history (setting only)
  const searchHistoryContext = React.useContext(SearchHistoryContext);
  const { updateSearchHistory } = searchHistoryContext || {};

  React.useEffect(() => {
    if (drugInfo.length > 0 && updateSearchHistory) {
      updateSearchHistory(params);
    }
  }, [drugInfo, params, updateSearchHistory]);

  const [showNavbar, setShowNavbar] = React.useState(false);
  const searchBoxSectionRef = React.useRef<HTMLDivElement>(null);
  const isSearchBoxInView = useInView(searchBoxSectionRef);

  React.useEffect(() => {
    setShowNavbar(!!!isSearchBoxInView);
  }, [isSearchBoxInView]);

  return (
    <>
      <AnimatePresence>
        {showNavbar && <Scrollbar key={"scrollbar"} />}
      </AnimatePresence>

      <InfoCard />

      <div
        ref={searchBoxSectionRef}
        className="d-flex flex-column justify-content-center align-items-center"
      >
        <Row className={"mb-4 mt-5"}>
          <Header />
        </Row>
        <Row className={"mb-4 text-center"}>
          <SearchBox />
        </Row>
      </div>

      {loading ? (
        <PropagateLoader
          className={"d-flex justify-content-center mt-4"}
          color={primaryColor}
          loading={loading}
          size={15}
        />
      ) : (
        <>
          {drugInfoError ? (
            <Col
              className={
                "d-flex flex-column justify-content-center align-items-center"
              }
            >
              <Bug
                weight={"light"}
                className={"display-1 text-secondary my-3"}
              />
              <h1 className={"display-1 mb-3"}>Oops!</h1>
              <h2 className={"text-secondary"}>
                We couldn't find what you were looking for...
              </h2>
            </Col>
          ) : (
            <Row className={"justify-content-center mx-auto"}>
              {params.searchBy === "brand_name"
                ? Object.entries(groupedByBrandName).map(
                    ([brandName, drugs], index) => (
                      <React.Fragment key={index}>
                        <Col
                          xs={12}
                          className={
                            "mb-1 d-flex flex-column justify-content-center align-items-center"
                          }
                        >
                          <div>
                            <h1>{brandName}</h1>
                            <hr />
                          </div>
                        </Col>
                        {(drugs as DrugProperties[]).map((drug, index) => (
                          <Col
                            xs={isMobile ? 12 : drugInfo.length === 1 ? 6 : 4}
                            key={index}
                            className="mb-4"
                          >
                            <DrugPropertyBox
                              drug={drug}
                              isSingle={drugs.length === 1}
                            />
                          </Col>
                        ))}
                      </React.Fragment>
                    ),
                  )
                : drugInfo.map((drug, index) => (
                    <Col
                      xs={isMobile ? 12 : drugInfo.length === 1 ? 6 : 4}
                      key={index}
                      className="mb-4"
                    >
                      <DrugPropertyBox
                        drug={drug}
                        isSingle={drugInfo.length === 1}
                      />
                    </Col>
                  ))}
            </Row>
          )}

          <Row className={"justify-content-center mx-auto"}>
            <ChartSection />
          </Row>

          <Row className={"justify-content-center mx-auto"}>
            <Col xs={isMobile ? 12 : 8}>
              <RelevantArticles />
            </Col>
          </Row>
        </>
      )}
    </>
  );
};

export default Result;
